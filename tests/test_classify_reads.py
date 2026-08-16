import pickle

import pandas as pd

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
