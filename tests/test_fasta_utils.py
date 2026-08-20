
from fasta_utils import parse_fasta
from fasta_utils import parse_fastq_qualities
from fasta_utils import extract_kmers
from fasta_utils import reverse_complement
from fasta_utils import canonical_kmer
from fasta_utils import extract_canonical_kmers
from classifier_functions import build_kmer_index

def test_parse_fasta_multiple_sequences(tmp_path):
    fake_fasta_content = ">seq1\nACGT\nACGT\n>seq2\nTTTT\n"
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fake_fasta_content)

    result = parse_fasta(str(fasta_file))

    assert result == {"seq1": "ACGTACGT", "seq2": "TTTT"}

def test_parse_fasta_single_sequence(tmp_path):
    fake_fasta_content = ">seq1\nACGT\nACGT\n"
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fake_fasta_content)

    result = parse_fasta(str(fasta_file))

    assert result == {"seq1": "ACGTACGT"}

def test_parse_fasta_empty_file(tmp_path):
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text("")

    result = parse_fasta(str(fasta_file))

    assert result == {}

def test_extract_kmers():
    sequence = "ACGTACGT"
    k = 3
    expected_kmers = ["ACG", "CGT", "GTA", "TAC", "ACG", "CGT"]
    
    result = extract_kmers(sequence, k)
    
    assert result == expected_kmers

def test_extract_kmers_k_greater_than_sequence_length():
    sequence = "ACGT"
    k = 5
    
    result = extract_kmers(sequence, k)
    
    assert result == []

def test_extract_kmers_k_is_zero():
    sequence = "ACGT"
    k = 0
    
    try:
        extract_kmers(sequence, k)
    except ValueError as e:
        assert str(e) == "k must be a positive integer"
    else:
        assert False, "Expected ValueError for k=0"

def test_extract_kmers_k_is_negative():
    sequence = "ACGT"
    k = -1
    
    try:
        extract_kmers(sequence, k)
    except ValueError as e:
        assert str(e) == "k must be a positive integer"
    else:
        assert False, "Expected ValueError for negative k"

def test_extract_kmers_k_equals_sequence_length():

    sequence = "ACGT"
    k = 4
    expected_kmers = ["ACGT"]
    
    result = extract_kmers(sequence, k)
    
    assert result == expected_kmers

def test_extract_kmers_k_equals_one():

    sequence = "ACGT"
    k = 1
    expected_kmers = ["A", "C", "G", "T"]
    
    result = extract_kmers(sequence, k)
    
    assert result == expected_kmers

# --- Reverse complement and canonical k-mers -------------------------------
#
# DNA is double-stranded and a sequencer reports whichever strand it happened
# to read, without saying which. The same physical fragment can therefore
# arrive written either way round. Canonicalizing every k-mer to
# min(kmer, revcomp(kmer)) collapses both forms onto one dictionary key, which
# is what lets the classifier match a read regardless of strand.
#
# Measured against the pre-canonical forward-only index on SRR10391187: a read
# and its reverse complement got DIFFERENT species for 0.66% of R1 reads and
# 5.84% of R2 reads, concentrated almost entirely in E. coli vs S. enterica.


def test_reverse_complement_basic():
    # Complement each base, then reverse the order (sequences are written 5'->3').
    assert reverse_complement("AAGC") == "GCTT"


def test_reverse_complement_is_its_own_inverse():
    """Applying it twice must return the original - it is a round trip."""
    sequence = "ACGTTGCAAGGCTTACCA"

    assert reverse_complement(reverse_complement(sequence)) == sequence


def test_reverse_complement_leaves_ambiguous_bases_alone():
    """N means 'unknown base', and the complement of unknown is unknown.

    ~2.3% of SRR10391187 reads contain at least one N, so this is not a
    hypothetical: crashing or dropping them would lose real reads.
    """
    assert reverse_complement("ACGTN") == "NACGT"


def test_reverse_complement_handles_lowercase():
    """Soft-masked FASTA (lowercase = repeat region) is common in public data."""
    assert reverse_complement("acgt") == "acgt"
    assert reverse_complement("aagc") == "gctt"


def test_reverse_complement_empty_sequence():
    assert reverse_complement("") == ""


def test_canonical_kmer_agrees_across_strands():
    """The whole point: both strands of one fragment reduce to one key."""
    kmer = "CGTT"
    rc = reverse_complement(kmer)

    assert rc != kmer, "test setup needs a k-mer that is not self-complementary"
    assert canonical_kmer(kmer) == canonical_kmer(rc)


def test_canonical_kmer_returns_one_of_the_two_forms():
    """It picks a representative, it does not invent a new string."""
    kmer = "CGTT"

    assert canonical_kmer(kmer) in (kmer, reverse_complement(kmer))


def test_canonical_kmer_is_idempotent():
    """Canonicalizing an already-canonical k-mer must be a no-op.

    If this failed, applying the transform twice anywhere in the pipeline
    would silently change keys.
    """
    kmer = "CGTT"
    once = canonical_kmer(kmer)

    assert canonical_kmer(once) == once


def test_canonical_kmer_self_complementary_input():
    """At even k a k-mer can equal its own reverse complement (ACGT does).

    Both branches of the comparison give the same answer, so it returns
    itself. At the real k=21 this case cannot arise: the middle base would
    have to pair with itself.
    """
    assert reverse_complement("ACGT") == "ACGT"
    assert canonical_kmer("ACGT") == "ACGT"


def test_extract_canonical_kmers_matches_manual_canonicalization():
    sequence = "ACGTTGCA"
    k = 4

    result = extract_canonical_kmers(sequence, k)

    assert result == [canonical_kmer(kmer) for kmer in extract_kmers(sequence, k)]
    assert result == ["ACGT", "AACG", "CAAC", "GCAA", "TGCA"]


def test_extract_canonical_kmers_preserves_position_count():
    """Canonicalizing rewrites k-mers; it never adds or drops positions.

    classify_read_top_hit's confidence denominator counts positions with
    plain extract_kmers, which is only valid because these two agree.
    """
    sequence = "ACGTTGCAAGGCTTACCA"

    assert len(extract_canonical_kmers(sequence, 4)) == len(extract_kmers(sequence, 4))


def test_extract_canonical_kmers_is_strand_invariant():
    """THE invariant this whole change exists to establish.

    A sequence and its reverse complement are the same physical molecule, so
    they must produce the same multiset of keys. Compared as sorted lists
    rather than in order, because reverse-complementing also reverses the
    order of the positions. That is fine for classification, which only
    tallies k-mers and never uses their order.
    """
    sequence = "ACGTTGCAAGGCTTACCA"

    forward = sorted(extract_canonical_kmers(sequence, 5))
    reverse = sorted(extract_canonical_kmers(reverse_complement(sequence), 5))

    assert forward == reverse


def test_extract_canonical_kmers_inherits_k_too_large():
    """Delegating to extract_kmers means this case is handled without a guard."""
    assert extract_canonical_kmers("ACGT", 5) == []


def test_extract_canonical_kmers_inherits_k_validation():
    try:
        extract_canonical_kmers("ACGT", 0)
    except ValueError as e:
        assert str(e) == "k must be a positive integer"
    else:
        assert False, "Expected ValueError for k=0"


def test_parse_fastq_qualities_maps_each_read(tmp_path):
    fake_fastq_content = "@read1\nACGT\n+\nIIII\n@read2\nTTTT\n+\n!!!!\n"
    fastq_file = tmp_path / "test.fastq"
    fastq_file.write_text(fake_fastq_content)

    result = parse_fastq_qualities(str(fastq_file))

    assert result == {"read1": "IIII", "read2": "!!!!"}


def test_parse_fastq_qualities_empty_file(tmp_path):
    fastq_file = tmp_path / "test.fastq"
    fastq_file.write_text("")

    result = parse_fastq_qualities(str(fastq_file))

    assert result == {}


def test_build_kmer_index_deduplicates():
    sequences = {
        "seq1": "ACGTACGT",
        "seq2": "TTTT"
    }
    k = 3
    expected_index = {
        "ACG": {"seq1"},
        "CGT": {"seq1"},
        "GTA": {"seq1"},
        "TAC": {"seq1"},
        "TTT": {"seq2"}
    }

    result = build_kmer_index(sequences, k)

    assert result == expected_index
