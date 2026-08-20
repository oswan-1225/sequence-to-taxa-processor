import pandas as pd
import pytest

from classifier_functions import build_kmer_index, save_index
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
    save_index(str(index_path), index, k)
    return str(index_path)


# Matches neither species in make_index(), on either strand. Since the index
# became canonical, "matches nothing" is a stronger condition than it looks:
# a read also has to miss the reverse complement of every reference k-mer, so
# a hand-picked run of one base will not do it (see the poly-C test below).
NO_MATCH_READ = "TTAGTTGTGCCG"


def test_classify_file_writes_expected_csv(tmp_path):
    index_path = make_index(tmp_path)
    reads_path = tmp_path / "reads.fasta"
    write_fasta(reads_path, [
        ("read1", "ACGTACGTACGT"),  # matches species_a
        ("read2", "TTTTGGGGTTTT"),  # matches species_b
        ("read3", NO_MATCH_READ),   # matches nothing, on either strand
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


def test_classify_file_matches_read_on_the_opposite_strand(tmp_path):
    """A read written on the reverse strand must reach the same species.

    species_b's genome contains "GGGG". A poly-C read is the reverse
    complement of a poly-G stretch, i.e. literally the other strand of the
    same duplex, so it belongs to species_b. Against the old forward-only
    index this read was reported as matching nothing at all.

    This is the unit-scale version of what was measured on SRR10391187,
    where a read and its reverse complement were assigned to different
    species for 5.84% of R2 reads.
    """
    index_path = make_index(tmp_path)
    reads_path = tmp_path / "reads.fasta"
    write_fasta(reads_path, [
        ("forward", "GGGGGGGGGGGG"),
        ("reverse", "CCCCCCCCCCCC"),
    ])
    output_path = tmp_path / "results.csv"

    results_df = classify_file(index_path, str(reads_path), k=4, output_path=str(output_path))

    by_read = results_df.set_index("read_id")
    assert by_read.loc["forward", "best_match"] == "species_b"
    assert by_read.loc["reverse", "best_match"] == "species_b"
    assert by_read.loc["forward", "confidence"] == by_read.loc["reverse", "confidence"]


# --- n_species_hit ----------------------------------------------------------
#
# best_match and confidence alone collapse four genuinely different outcomes
# into two indistinguishable pairs on disk: an unambiguous single-species hit
# looks identical to a winner-take-all one, and a tie looks identical to a read
# that matched nothing at all.
#
# That cost two reported numbers their reproducibility. "Discard ambiguous
# reads" accuracy (README: 138.7%) needs to select reads that hit exactly one
# species, and the tie rate (README: 0.15%) needs to separate ties from
# no-matches. Both were measured by throwaway scripts and neither could be
# recomputed from any saved run. n_species_hit is the smallest thing that
# answers both.


def test_classify_file_records_all_four_ambiguity_cases(tmp_path):
    """n_species_hit + best_match must distinguish all four outcomes."""
    index_path = make_index(tmp_path)
    reads_path = tmp_path / "reads.fasta"
    write_fasta(reads_path, [
        ("unambiguous", "ACGTACGTACGT"),    # hits species_a only
        ("ambiguous", "ACGTGGGGTTTT"),      # hits both, species_b wins outright
        ("tied", "ACGTACGTTTTGGG"),         # hits both, 5 votes each
        ("no_match", NO_MATCH_READ),        # hits nothing, on either strand
    ])
    output_path = tmp_path / "results.csv"

    results_df = classify_file(index_path, str(reads_path), k=4, output_path=str(output_path))
    by_read = results_df.set_index("read_id")

    # Unambiguous: exactly one species, and it won.
    assert by_read.loc["unambiguous", "n_species_hit"] == 1
    assert by_read.loc["unambiguous", "best_match"] == "species_a"

    # Ambiguous: several species in play, one took the whole read.
    assert by_read.loc["ambiguous", "n_species_hit"] == 2
    assert by_read.loc["ambiguous", "best_match"] == "species_b"

    # Tie: several species in play, none won. Distinguished from no_match
    # ONLY by n_species_hit - both have a null best_match.
    assert by_read.loc["tied", "n_species_hit"] == 2
    assert pd.isna(by_read.loc["tied", "best_match"])

    # No match: nothing in the reference set resembles this read.
    assert by_read.loc["no_match", "n_species_hit"] == 0
    assert pd.isna(by_read.loc["no_match", "best_match"])


def test_n_species_hit_makes_discarded_and_tie_rates_recomputable(tmp_path):
    """The two README numbers, derived from a saved file rather than a script."""
    index_path = make_index(tmp_path)
    reads_path = tmp_path / "reads.fasta"
    write_fasta(reads_path, [
        ("unambiguous", "ACGTACGTACGT"),
        ("ambiguous", "ACGTGGGGTTTT"),
        ("tied", "ACGTACGTTTTGGG"),
        ("no_match", NO_MATCH_READ),
    ])
    output_path = tmp_path / "results.csv"
    classify_file(index_path, str(reads_path), k=4, output_path=str(output_path))

    on_disk = pd.read_csv(output_path)

    # "Discard ambiguous multi-genome reads" keeps only single-species hits.
    unambiguous = on_disk[on_disk["n_species_hit"] == 1]
    assert list(unambiguous["read_id"]) == ["unambiguous"]

    # A tie is a null best_match that still hit something.
    ties = on_disk[on_disk["best_match"].isna() & (on_disk["n_species_hit"] >= 2)]
    assert list(ties["read_id"]) == ["tied"]

    # A no-match is a null best_match that hit nothing.
    no_match = on_disk[on_disk["best_match"].isna() & (on_disk["n_species_hit"] == 0)]
    assert list(no_match["read_id"]) == ["no_match"]


def test_n_species_hit_survives_the_csv_round_trip(tmp_path):
    index_path = make_index(tmp_path)
    reads_path = tmp_path / "reads.fasta"
    write_fasta(reads_path, [("read1", "ACGTGGGGTTTT")])
    output_path = tmp_path / "results.csv"

    results_df = classify_file(index_path, str(reads_path), k=4, output_path=str(output_path))
    on_disk = pd.read_csv(output_path)

    assert "n_species_hit" in on_disk.columns
    assert on_disk.loc[0, "n_species_hit"] == results_df.loc[0, "n_species_hit"]


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
        ("read3", NO_MATCH_READ),   # matches nothing, on either strand
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
