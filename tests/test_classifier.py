from classifier import classify_read_top_hit
from classifier import build_kmer_index

def test_classify_read_top_hit_clear_match():
    index = build_kmer_index({"species_a": "ACGTACGT", "species_b": "TTTTGGGG"}, k=4)
    result = classify_read_top_hit("ACGTACGT", index, k=4)

    assert result['best_match'] == "species_a"
    assert result['confidence'] > 0.5

def test_classify_read_top_hit_no_match():
    index = build_kmer_index({"species_a": "ACGTACGT"}, k=4)
    result = classify_read_top_hit("TTTTTTTT", index, k=4)

    assert result['best_match'] is None
    assert result['confidence'] == 0.0

def test_classify_read_top_hit_read_shorter_than_k():
    index = build_kmer_index({"species_a": "ACGTACGT"}, k=4)
    result = classify_read_top_hit("AC", index, k=4)

    assert result['best_match'] is None
    assert result['confidence'] == 0.0