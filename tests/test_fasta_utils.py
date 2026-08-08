
from fasta_utils import parse_fasta

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