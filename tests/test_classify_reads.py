import pickle

import pandas as pd
import pytest

from classifier_functions import build_kmer_index
from classify_reads import classify_file
from database import get_abundance


def write_fasta(path, records):
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n{seq}\n")


def make_index(tmp_path, k=4):
    genomes = {
        "species_a": "ACGTACGTACGT",
        "species_b": "TTTTGGGGTTTT",
    }
    index = build_kmer_index(genomes, k)
    index_path = tmp_path / "index.pkl"
    with open(index_path, "wb") as f:
        pickle.dump(index, f)
    return str(index_path)


def test_classify_file_writes_expected_csv(tmp_path):
    index_path = make_index(tmp_path)
    reads_path = tmp_path / "reads.fasta"
    write_fasta(reads_path, [
        ("read1", "ACGTACGTACGT"),  # matches species_a
        ("read2", "TTTTGGGGTTTT"),  # matches species_b
        ("read3", "CCCCCCCCCCCC"),  # matches nothing
    ])
    output_path = tmp_path / "results.csv"

    results_df = classify_file(index_path, str(reads_path), k=4, output_path=str(output_path))

    assert output_path.exists()
    on_disk = pd.read_csv(output_path)
    expected = results_df.reset_index(drop=True).fillna(pd.NA)
    pd.testing.assert_frame_equal(expected, on_disk.fillna(pd.NA))

    by_read = results_df.set_index("read_id")
    assert by_read.loc["read1", "best_match"] == "species_a"
    assert by_read.loc["read2", "best_match"] == "species_b"
    assert pd.isna(by_read.loc["read3", "best_match"])
    assert by_read.loc["read3", "confidence"] == 0.0


def test_classify_file_min_quality_drops_low_quality_reads(tmp_path):
    index_path = make_index(tmp_path)
    reads_path = tmp_path / "reads.fastq"
    with open(reads_path, "w") as f:
        f.write("@read1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n")  # mean Q40 - kept
        f.write("@read2\nTTTTGGGGTTTT\n+\n!!!!!!!!!!!!\n")  # mean Q0 - dropped
    output_path = tmp_path / "results.csv"

    results_df = classify_file(index_path, str(reads_path), k=4, output_path=str(output_path),
                                min_quality=20.0)

    assert set(results_df["read_id"]) == {"read1"}


def test_classify_file_min_quality_rejects_fasta_input(tmp_path):
    index_path = make_index(tmp_path)
    reads_path = tmp_path / "reads.fasta"
    write_fasta(reads_path, [("read1", "ACGTACGTACGT")])
    output_path = tmp_path / "results.csv"

    try:
        classify_file(index_path, str(reads_path), k=4, output_path=str(output_path), min_quality=20.0)
    except ValueError as e:
        assert "FASTQ" in str(e)
    else:
        assert False, "Expected ValueError for min_quality on FASTA input"


def test_classify_file_redistribute_off_by_default(tmp_path):
    index_path = make_index(tmp_path)
    reads_path = tmp_path / "reads.fasta"
    write_fasta(reads_path, [("read1", "ACGTACGTACGT")])
    output_path = tmp_path / "results.csv"

    classify_file(index_path, str(reads_path), k=4, output_path=str(output_path))

    redistributed_path = tmp_path / "results_redistributed.csv"
    assert not redistributed_path.exists()


def test_classify_file_redistribute_splits_ambiguous_votes(tmp_path):
    index_path = make_index(tmp_path)
    reads_path = tmp_path / "reads.fasta"
    write_fasta(reads_path, [
        ("read1", "ACGTACGTACGT"),  # pure species_a match
        ("read2", "ACGTGGGGTTTT"),  # ambiguous mix of both
        ("read3", "CCCCCCCCCCCC"),  # matches nothing
    ])
    output_path = tmp_path / "results.csv"

    classify_file(index_path, str(reads_path), k=4, output_path=str(output_path), redistribute=True)

    redistributed_path = tmp_path / "results_redistributed.csv"
    assert redistributed_path.exists()
    redistributed_df = pd.read_csv(redistributed_path).set_index("species")
    # 2 classified reads total: read1 gives species_a full credit, read2
    # splits its credit across both species it hit.
    assert redistributed_df["estimated_reads"].sum() == pytest.approx(2.0)
    assert redistributed_df.loc["species_a", "estimated_reads"] > redistributed_df.loc["species_b", "estimated_reads"]


def test_classify_file_populates_database(tmp_path):
    index_path = make_index(tmp_path)
    reads_path = tmp_path / "SRR000001_1.fastq"
    with open(reads_path, "w") as f:
        f.write("@read1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n")
    output_path = tmp_path / "results.csv"
    db_path = tmp_path / "results.db"

    classify_file(index_path, str(reads_path), k=4, output_path=str(output_path),
                   db_path=str(db_path), source="test")

    abundance_df = get_abundance(str(db_path))
    row = abundance_df.iloc[0]
    assert row["sample_id"] == "SRR000001"
    assert row["best_match"] == "species_a"
    assert row["count"] == 1
