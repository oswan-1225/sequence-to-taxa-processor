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


# --- n_species_hit persistence and schema versioning ------------------------


def test_insert_sample_results_persists_n_species_hit(tmp_path):
    db_path = str(tmp_path / "test.db")
    create_database(db_path)
    results = [
        {"read_id": "r1", "best_match": "species_a", "confidence": 1.0, "n_species_hit": 1},
        {"read_id": "r2", "best_match": "species_a", "confidence": 0.6, "n_species_hit": 3},
        {"read_id": "r3", "best_match": None, "confidence": 0.0, "n_species_hit": 2},
        {"read_id": "r4", "best_match": None, "confidence": 0.0, "n_species_hit": 0},
    ]
    insert_sample_results(db_path, "sample1", "local", results)

    conn = sqlite3.connect(db_path)
    rows = dict(conn.execute(
        "SELECT read_id, n_species_hit FROM classifications"
    ).fetchall())
    conn.close()

    assert rows == {"r1": 1, "r2": 3, "r3": 2, "r4": 0}


def test_ties_and_no_matches_are_distinguishable_in_the_database(tmp_path):
    """Both have a NULL best_match; only n_species_hit separates them."""
    db_path = str(tmp_path / "test.db")
    create_database(db_path)
    insert_sample_results(db_path, "sample1", "local", [
        {"read_id": "tied", "best_match": None, "confidence": 0.0, "n_species_hit": 2},
        {"read_id": "nothing", "best_match": None, "confidence": 0.0, "n_species_hit": 0},
    ])

    conn = sqlite3.connect(db_path)
    ties = conn.execute(
        "SELECT read_id FROM classifications "
        "WHERE best_match IS NULL AND n_species_hit >= 2"
    ).fetchall()
    conn.close()

    assert ties == [("tied",)]


def test_insert_tolerates_results_without_n_species_hit(tmp_path):
    """Read with .get(), so hand-built result dicts predating the column work."""
    db_path = str(tmp_path / "test.db")
    create_database(db_path)

    insert_sample_results(db_path, "sample1", "local", [
        {"read_id": "r1", "best_match": "species_a", "confidence": 1.0},
    ])

    conn = sqlite3.connect(db_path)
    value = conn.execute("SELECT n_species_hit FROM classifications").fetchone()[0]
    conn.close()

    assert value is None


def test_create_database_rejects_outdated_schema(tmp_path):
    """CREATE TABLE IF NOT EXISTS is a no-op on an existing table.

    Without an explicit check, an old database keeps its old columns and the
    failure surfaces much later as an opaque sqlite error inside executemany.
    """
    db_path = str(tmp_path / "old.db")
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE classifications (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            sample_id TEXT,
            read_id TEXT,
            best_match TEXT,
            confidence REAL
        )
    """)
    conn.commit()
    conn.close()

    with pytest.raises(ValueError, match="outdated classifications table"):
        create_database(db_path)


def test_outdated_schema_error_names_the_missing_column_and_the_fix(tmp_path):
    db_path = str(tmp_path / "old.db")
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE classifications (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            sample_id TEXT,
            read_id TEXT,
            best_match TEXT,
            confidence REAL
        )
    """)
    conn.commit()
    conn.close()

    with pytest.raises(ValueError) as excinfo:
        create_database(db_path)

    message = str(excinfo.value)
    assert "n_species_hit" in message
    assert "delete it and rerun" in message


def test_create_database_is_still_idempotent_on_a_current_database(tmp_path):
    """The guard must not fire on a database this version wrote."""
    db_path = str(tmp_path / "test.db")
    create_database(db_path)

    create_database(db_path)  # must not raise
