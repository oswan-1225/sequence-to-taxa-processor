import itertools
import os
import subprocess
import sys

from classifier_functions import classify_read_top_hit
from classifier_functions import build_kmer_index
from fasta_utils import extract_canonical_kmers
from fasta_utils import reverse_complement

SRC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src")

def test_classify_read_top_hit_clear_match():
    index = build_kmer_index({"species_a": "ACGTACGT", "species_b": "TTTTGGGG"}, k=4)
    result = classify_read_top_hit("ACGTACGT", index, k=4)

    assert result['best_match'] == "species_a"
    assert result['confidence'] > 0.5
    assert result['votes'] == {"species_a": 5}

def test_classify_read_top_hit_no_match():
    index = build_kmer_index({"species_a": "ACGTACGT"}, k=4)
    result = classify_read_top_hit("TTTTTTTT", index, k=4)

    assert result['best_match'] is None
    assert result['confidence'] == 0.0
    assert result['votes'] == {}

def test_classify_read_top_hit_read_shorter_than_k():
    index = build_kmer_index({"species_a": "ACGTACGT"}, k=4)
    result = classify_read_top_hit("AC", index, k=4)

    assert result['best_match'] is None
    assert result['confidence'] == 0.0
    assert result['votes'] == {}

def test_classify_read_top_hit_ambiguous_match():
    index = build_kmer_index({"species_a": "ACGTACGT", "species_b": "ACGTTTTT"}, k=4)
    result = classify_read_top_hit("ACGTACGT", index, k=4)

    assert result['best_match'] == "species_a"
    # Both species score matched k-mers here, but confidence alone (1.0)
    # would hide that species_b was a real competitor - only 'votes' shows
    # species_b got credit too. This is the exact gap the abundance
    # re-estimation work needs 'votes' to close.
    assert result['votes'] == {"species_a": 5, "species_b": 2}


# --- Tie handling and reproducibility -------------------------------------
#
# build_kmer_index stores each k-mer's species in a set, and set iteration
# order for strings varies between processes because Python randomizes string
# hashing. A tied read therefore used to be awarded to whichever species that
# randomized order happened to put first, so the same input produced different
# output on different runs. Ties now resolve to best_match None instead.
#
# The trap these tests are built around: set iteration order is FIXED within a
# single process. Calling the classifier twice in one test proves nothing, and
# that is why the pre-existing suite never caught this. Each test below breaks
# out of the single-process assumption in a different way.


# Every CANONICAL 4-mer in this read is distinct. Tests that give one species
# a hand-tuned vote margin need that: in a read like "ACGTACGTACGT" the k-mer
# "ACGT" recurs, so patching a single index entry would silently change
# several votes at once instead of one.
#
# Canonicalization makes this constraint stricter than it used to be. The
# previous constant here, "ACGTTGCAAGGCTTACCA", has 15 distinct forward 4-mers
# but only 14 distinct canonical ones, because two of its positions are
# reverse complements of each other and now collapse onto one key.
UNIQUE_KMER_READ = "ACCCCATCGGACTGGCAT"


def _index_with_order(read, k, species_order):
    """Build an index by hand where every k-mer maps to species in a chosen order.

    classify_read only iterates kmer_index[kmer], so a list works in place of a
    set and gives the test explicit control over the order votes are tallied in
    - the thing the real set makes unpredictable.

    Keys must be CANONICAL, because that is what classify_read looks up. Keying
    on raw forward k-mers does not fail loudly: the positions that happen to
    already be canonical still match, so the test keeps passing while quietly
    exercising only a subset of the read.
    """
    return {kmer: list(species_order) for kmer in extract_canonical_kmers(read, k)}


def test_tie_is_unclassified_regardless_of_vote_order():
    """Every permutation of a tied read's vote order must give the same answer."""
    read = "ACGTACGTACGTTTGACC"
    species = ["species_a", "species_b", "species_c"]

    for order in itertools.permutations(species):
        result = classify_read_top_hit(read, _index_with_order(read, 4, order), k=4)

        assert result['best_match'] is None, f"tie resolved to a winner for order {order}"
        assert result['confidence'] == 0.0
        # The tally survives even though no winner was picked - callers that
        # want to know what the read hit still can.
        assert set(result['votes']) == set(species)
        assert len(set(result['votes'].values())) == 1, "test setup is not actually a tie"


def test_clear_winner_is_unaffected_by_vote_order():
    """The 95% of reads with an outright winner must still classify normally.

    Guards the branch the tie fix did not intend to touch: an early version
    collected vote counts instead of species names, which left tie detection
    working perfectly while writing an integer into best_match.
    """
    read = UNIQUE_KMER_READ
    kmers = extract_canonical_kmers(read, 4)
    assert len(kmers) == len(set(kmers)), "read must have unique canonical k-mers for this setup"

    for order in (["species_a", "species_b"], ["species_b", "species_a"]):
        index = {kmer: ["species_a"] for kmer in kmers}
        index[kmers[0]] = list(order)  # species_b picks up exactly one vote

        result = classify_read_top_hit(read, index, k=4)

        assert result['best_match'] == "species_a"
        assert isinstance(result['best_match'], str)
        assert result['confidence'] == 1.0
        assert result['votes'] == {"species_a": len(kmers), "species_b": 1}


def test_single_vote_margin_still_picks_a_winner():
    """A one-vote lead is a winner, not a tie - the boundary the fix sits on."""
    read = UNIQUE_KMER_READ
    kmers = extract_canonical_kmers(read, 4)
    assert len(kmers) == len(set(kmers)), "read must have unique canonical k-mers for this setup"
    index = {kmer: ["species_a", "species_b"] for kmer in kmers}
    index[kmers[0]] = ["species_a"]  # a leads b by exactly one

    result = classify_read_top_hit(read, index, k=4)

    assert result['best_match'] == "species_a"
    assert result['votes']["species_a"] == result['votes']["species_b"] + 1


# --- Strand invariance ------------------------------------------------------
#
# DNA is double-stranded and a sequencer reports whichever strand it read,
# without saying which. A read and its reverse complement are the same physical
# molecule, so the classifier must return the same answer for both.
#
# The pre-canonical implementation did not. Measured on SRR10391187 against the
# forward-only index, a read and its reverse complement were assigned to
# DIFFERENT species for 0.66% of R1 reads and 5.84% of R2 reads, concentrated
# almost entirely in E. coli vs S. enterica - two Enterobacteriaceae sharing
# enough 16S that a strand flip tips the vote between them. It moved the
# community profile by ~1.5 percentage points on each of that pair.
#
# This now holds by construction rather than by approximation, since
# canonical_kmer(reverse_complement(x)) == canonical_kmer(x) for every k-mer.
# So these assert exact equality, not "close enough".


def test_read_and_its_reverse_complement_classify_identically():
    """The property the canonical-k-mer change exists to establish."""
    index = build_kmer_index({
        "species_a": "ACGTTGCAAGGCTTACCAGGTTACG",
        "species_b": "TTTTGGGGCCCCAAAATTTTGGGGC",
    }, k=5)
    read = "TTGCAAGGCTTACCAG"

    forward = classify_read_top_hit(read, index, k=5)
    reverse = classify_read_top_hit(reverse_complement(read), index, k=5)

    assert forward['best_match'] == reverse['best_match'] == "species_a"
    assert forward['votes'] == reverse['votes']
    assert forward['confidence'] == reverse['confidence']


def test_read_matches_reference_stored_on_the_opposite_strand():
    """A genome that carries a region in reverse orientation still matches.

    This is why the forward-only index did not simply leave reverse-strand
    reads unclassified: real genomes transcribe genes from both strands, so
    both orientations are usually present somewhere in the reference set. The
    failure mode was a wrong answer, not a missing one.
    """
    region = "ACGTTGCAAGGCTTACCAGG"
    index = build_kmer_index({
        "forward_species": "GGGGGG" + region + "GGGGGG",
        "reverse_species": "AAAAAA" + reverse_complement(region) + "AAAAAA",
    }, k=7)

    result = classify_read_top_hit(region, index, k=7)

    # The region is present in both genomes, one on each strand. A
    # strand-agnostic classifier sees a genuine tie rather than picking the
    # species whose orientation happens to match how the read was written.
    assert result['best_match'] is None
    assert set(result['votes']) == {"forward_species", "reverse_species"}
    assert len(set(result['votes'].values())) == 1


def test_strand_invariance_holds_for_unclassified_reads_too():
    """A read that matches nothing must match nothing on both strands."""
    index = build_kmer_index({"species_a": "ACGTACGTACGTACGT"}, k=5)
    read = "TTAGTTGTGCCGATAC"

    forward = classify_read_top_hit(read, index, k=5)
    reverse = classify_read_top_hit(reverse_complement(read), index, k=5)

    assert forward['best_match'] is None
    assert forward == reverse


def test_classification_is_reproducible_across_hash_seeds():
    """End-to-end guard: identical output under different string-hash seeds.

    This is the test that would have caught the original bug. It runs the real
    set-backed index in fresh interpreters, since PYTHONHASHSEED can only be
    set at process start - it cannot be changed from inside a running test.
    """
    program = (
        "import sys; sys.path.insert(0, %r)\n"
        "from classifier_functions import build_kmer_index, classify_read_top_hit\n"
        "read = 'ACGTACGTACGTTTGACC'\n"
        "genomes = {n: read for n in "
        "['Species_alpha', 'Species_beta', 'Species_gamma', 'Species_delta']}\n"
        "index = build_kmer_index(genomes, 4)\n"
        "print(classify_read_top_hit(read, index, 4)['best_match'])\n"
    ) % SRC_DIR

    outputs = set()
    for seed in ("0", "1", "2", "3", "4"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        completed = subprocess.run(
            [sys.executable, "-c", program],
            capture_output=True, text=True, env=env, check=True,
        )
        outputs.add(completed.stdout.strip())

    assert outputs == {"None"}, f"output varied with hash seed: {outputs}"