import os

import pandas as pd
import pytest

from classifier_functions import build_kmer_index, save_index
from database import get_classification_totals
from pipeline import run_pipeline


def write_fasta(path, records):
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n{seq}\n")


def make_genome_dir(tmp_path):
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    write_fasta(genome_dir / "species_a.fna", [("chr1", "ACGTACGTACGT")])
    write_fasta(genome_dir / "species_b.fna", [("chr1", "TTTTGGGGTTTT")])
    return str(genome_dir)


def make_reads_file(tmp_path, name="reads.fasta"):
    reads_path = tmp_path / name
    write_fasta(reads_path, [
        ("read1", "ACGTACGTACGT"),  # matches species_a
        ("read2", "TTTTGGGGTTTT"),  # matches species_b
    ])
    return str(reads_path)


def test_run_pipeline_produces_all_expected_outputs(tmp_path):
    genome_dir = make_genome_dir(tmp_path)
    reads_path = make_reads_file(tmp_path)
    output_dir = tmp_path / "out"

    outputs = run_pipeline(genome_dir, reads_path, str(output_dir), k=4)

    assert set(outputs.keys()) == {
        "index", "classifications_csv", "database", "diversity_csv", "species_abundance_csv", "plot"
    }
    for path in outputs.values():
        assert os.path.exists(path)


def test_run_pipeline_skip_plot_omits_plot_output(tmp_path):
    genome_dir = make_genome_dir(tmp_path)
    reads_path = make_reads_file(tmp_path)
    output_dir = tmp_path / "out"

    outputs = run_pipeline(genome_dir, reads_path, str(output_dir), k=4, skip_plot=True)

    assert "plot" not in outputs
    assert not (output_dir / "abundance_plot.png").exists()


def test_run_pipeline_reuses_existing_index_and_skips_build(tmp_path):
    reads_path = make_reads_file(tmp_path)
    output_dir = tmp_path / "out"

    genomes = {"species_a": "ACGTACGTACGT", "species_b": "TTTTGGGGTTTT"}
    index_path = tmp_path / "prebuilt.pkl"
    save_index(str(index_path), build_kmer_index(genomes, 4), 4)

    outputs = run_pipeline(None, reads_path, str(output_dir), k=4, index=str(index_path))

    assert outputs["index"] == str(index_path)
    assert not (output_dir / "kmer_index.pkl").exists()


def test_run_pipeline_requires_index_or_genome_dir(tmp_path):
    reads_path = make_reads_file(tmp_path)

    with pytest.raises(ValueError, match="Either --index or --genome-dir"):
        run_pipeline(None, reads_path, str(tmp_path / "out"), k=4)


def test_run_pipeline_rejects_missing_index_path(tmp_path):
    reads_path = make_reads_file(tmp_path)

    with pytest.raises(ValueError, match="does not exist"):
        run_pipeline(None, reads_path, str(tmp_path / "out"), k=4, index=str(tmp_path / "missing.pkl"))


def test_run_pipeline_fails_fast_and_stops_before_later_stages(tmp_path):
    genome_dir = make_genome_dir(tmp_path)
    output_dir = tmp_path / "out"
    missing_reads_path = str(tmp_path / "does_not_exist.fasta")

    with pytest.raises(RuntimeError, match="Stage 2 \\(classify reads\\) failed"):
        run_pipeline(genome_dir, missing_reads_path, str(output_dir), k=4)

    # the index was built (stage 1 succeeded), but nothing from stage 3/4 exists
    assert (output_dir / "kmer_index.pkl").exists()
    assert not (output_dir / "diversity_report.csv").exists()
    assert not (output_dir / "abundance_plot.png").exists()


def test_run_pipeline_shares_db_across_samples_without_clobbering(tmp_path):
    genome_dir = make_genome_dir(tmp_path)
    reads1 = make_reads_file(tmp_path, name="reads.fasta")
    reads2 = make_reads_file(tmp_path, name="reads2.fasta")
    shared_db = tmp_path / "shared.db"

    # shared_db doesn't exist yet - the first run must create it, not require
    # it to pre-exist.
    assert not shared_db.exists()
    outputs1 = run_pipeline(genome_dir, reads1, str(tmp_path / "out1"), k=4, db=str(shared_db))
    assert shared_db.exists()

    outputs2 = run_pipeline(genome_dir, reads2, str(tmp_path / "out2"), k=4, db=str(shared_db))

    assert outputs1["database"] == str(shared_db)
    assert outputs2["database"] == str(shared_db)

    # both samples' rows survive in the one shared database - the second
    # run must not have deleted the first sample's rows.
    all_sample_ids = set(get_classification_totals(str(shared_db))["sample_id"])
    assert all_sample_ids == {"reads", "reads2"}

    # each run's own report/plot outputs stay scoped to just that run's
    # sample, even though the database now holds two.
    abundance1 = pd.read_csv(outputs1["species_abundance_csv"])
    abundance2 = pd.read_csv(outputs2["species_abundance_csv"])
    assert set(abundance1["sample_id"]) == {"reads"}
    assert set(abundance2["sample_id"]) == {"reads2"}

    diversity1 = pd.read_csv(outputs1["diversity_csv"])
    diversity2 = pd.read_csv(outputs2["diversity_csv"])
    assert len(diversity1) == 1
    assert len(diversity2) == 1
    assert diversity1.iloc[0]["sample_id"] == "reads"
    assert diversity2.iloc[0]["sample_id"] == "reads2"


def test_run_pipeline_rerunning_same_sample_replaces_not_doubles(tmp_path):
    genome_dir = make_genome_dir(tmp_path)
    reads_path = make_reads_file(tmp_path)
    shared_db = tmp_path / "shared.db"

    run_pipeline(genome_dir, reads_path, str(tmp_path / "out1"), k=4, db=str(shared_db))
    run_pipeline(genome_dir, reads_path, str(tmp_path / "out2"), k=4, db=str(shared_db))

    totals = get_classification_totals(str(shared_db))
    assert len(totals) == 1
    assert totals.iloc[0]["sample_id"] == "reads"
