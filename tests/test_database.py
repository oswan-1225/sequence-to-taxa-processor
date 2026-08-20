import sqlite3

import pytest

from database import create_database, insert_sample_results, get_abundance, get_classification_totals


def test_create_database_creates_tables(tmp_path):
    db_path = str(tmp_path / "test.db")
    create_database(db_path)

    conn = sqlite3.connect(db_path)
    tables = {row[0] for row in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'"
    )}
    conn.close()

    assert {"samples", "classifications"} <= tables


def test_create_database_creates_missing_parent_dir(tmp_path):
    db_path = str(tmp_path / "nested" / "dir" / "test.db")
    create_database(db_path)

    assert (tmp_path / "nested" / "dir" / "test.db").exists()


def test_create_database_accepts_a_bare_filename(tmp_path, monkeypatch):
    """A db_path with no directory part must not crash.

    os.path.dirname("results.db") is "", and os.makedirs("") raises
    FileNotFoundError. Reachable via `classify_reads.py --db results.db`.
    Every other test passes an absolute tmp_path, so nothing else covers it.
    """
    monkeypatch.chdir(tmp_path)
    create_database("results.db")

    assert (tmp_path / "results.db").exists()


def test_insert_sample_results_inserts_sample_row_once(tmp_path):
    db_path = str(tmp_path / "test.db")
    create_database(db_path)
    results = [{"read_id": "r1", "best_match": "species_a", "confidence": 1.0}]

    insert_sample_results(db_path, "sample1", "local", results)
    insert_sample_results(db_path, "sample1", "local", results)

    conn = sqlite3.connect(db_path)
    count = conn.execute(
        "SELECT COUNT(*) FROM samples WHERE sample_id = ?", ("sample1",)
    ).fetchone()[0]
    classification_count = conn.execute(
        "SELECT COUNT(*) FROM classifications WHERE sample_id = ?", ("sample1",)
    ).fetchone()[0]
    conn.close()

    assert count == 1
    # The samples row was already deduplicated, but the classifications rows
    # were not: this second assertion is the one that catches a re-run
    # appending a duplicate set of reads.
    assert classification_count == 1


def test_reinserting_a_sample_replaces_its_classifications(tmp_path):
    """Re-running a sample must replace its rows, not append to them.

    pipeline.py always writes to output_dir/classifications.db, so running it
    twice into the same --output-dir used to double every abundance count
    while leaving the percentages unchanged, which is why nothing looked wrong.
    """
    db_path = str(tmp_path / "test.db")
    create_database(db_path)
    first = [{"read_id": f"r{i}", "best_match": "species_a", "confidence": 0.9}
             for i in range(5)]
    second = [{"read_id": f"r{i}", "best_match": "species_b", "confidence": 0.9}
              for i in range(3)]

    insert_sample_results(db_path, "sample1", "local", first)
    insert_sample_results(db_path, "sample1", "local", second)

    abundance = get_abundance(db_path)

    assert list(abundance["best_match"]) == ["species_b"], "old rows survived the re-run"
    assert int(abundance["count"].iloc[0]) == 3
    assert int(abundance["total"].iloc[0]) == 3


def test_reinserting_a_sample_leaves_other_samples_alone(tmp_path):
    """Replacement must be scoped per sample - the --db multi-sample flag needs this."""
    db_path = str(tmp_path / "test.db")
    create_database(db_path)
    sample1 = [{"read_id": "r1", "best_match": "species_a", "confidence": 0.9}]
    sample2 = [{"read_id": "r1", "best_match": "species_b", "confidence": 0.9},
               {"read_id": "r2", "best_match": "species_b", "confidence": 0.9}]

    insert_sample_results(db_path, "sample1", "local", sample1)
    insert_sample_results(db_path, "sample2", "local", sample2)
    insert_sample_results(db_path, "sample2", "local", sample2)  # re-run sample2 only

    totals = get_classification_totals(db_path).set_index("sample_id")

    assert int(totals.loc["sample1", "total_reads"]) == 1, "sample1 was collateral damage"
    assert int(totals.loc["sample2", "total_reads"]) == 2


def test_insert_sample_results_inserts_one_row_per_result(tmp_path):
    db_path = str(tmp_path / "test.db")
    create_database(db_path)
    results = [
        {"read_id": "r1", "best_match": "species_a", "confidence": 1.0},
        {"read_id": "r2", "best_match": "species_b", "confidence": 0.5},
        {"read_id": "r3", "best_match": None, "confidence": 0.0},
    ]

    insert_sample_results(db_path, "sample1", "local", results)

    conn = sqlite3.connect(db_path)
    rows = conn.execute(
        "SELECT read_id, best_match, confidence FROM classifications WHERE sample_id = ?",
        ("sample1",),
    ).fetchall()
    conn.close()

    assert set(rows) == {
        ("r1", "species_a", 1.0),
        ("r2", "species_b", 0.5),
        ("r3", None, 0.0),
    }


def test_get_abundance_counts_and_percentages(tmp_path):
    db_path = str(tmp_path / "test.db")
    create_database(db_path)
    results = [
        {"read_id": "r1", "best_match": "species_a", "confidence": 1.0},
        {"read_id": "r2", "best_match": "species_a", "confidence": 1.0},
        {"read_id": "r3", "best_match": "species_b", "confidence": 1.0},
        {"read_id": "r4", "best_match": None, "confidence": 0.0},
    ]
    insert_sample_results(db_path, "sample1", "local", results)

    abundance = get_abundance(db_path)

    row_a = abundance[abundance["best_match"] == "species_a"].iloc[0]
    row_b = abundance[abundance["best_match"] == "species_b"].iloc[0]

    assert row_a["count"] == 2
    assert row_a["total"] == 3  # unclassified read excluded from total
    assert row_a["percent"] == pytest.approx(200 / 3)
    assert row_b["count"] == 1
    assert row_b["percent"] == pytest.approx(100 / 3)


def test_get_classification_totals_counts_unclassified(tmp_path):
    db_path = str(tmp_path / "test.db")
    create_database(db_path)
    results = [
        {"read_id": "r1", "best_match": "species_a", "confidence": 1.0},
        {"read_id": "r2", "best_match": "species_a", "confidence": 1.0},
        {"read_id": "r3", "best_match": "species_b", "confidence": 1.0},
        {"read_id": "r4", "best_match": None, "confidence": 0.0},
    ]
    insert_sample_results(db_path, "sample1", "local", results)

    totals = get_classification_totals(db_path)
    row = totals.iloc[0]

    assert row["sample_id"] == "sample1"
    assert row["total_reads"] == 4
    assert row["classified_reads"] == 3
    assert row["unclassified_reads"] == 1
    assert row["unclassified_percent"] == pytest.approx(25.0)


def test_get_classification_totals_all_unclassified(tmp_path):
    db_path = str(tmp_path / "test.db")
    create_database(db_path)
    results = [{"read_id": "r1", "best_match": None, "confidence": 0.0}]
    insert_sample_results(db_path, "sample1", "local", results)

    totals = get_classification_totals(db_path)
    row = totals.iloc[0]

    assert row["classified_reads"] == 0
    assert row["unclassified_percent"] == pytest.approx(100.0)
