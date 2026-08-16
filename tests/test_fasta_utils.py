
from fasta_utils import parse_fasta
from fasta_utils import parse_fastq_qualities
from fasta_utils import extract_kmers
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
