from classifier_functions import classify_read_top_hit
from classifier_functions import build_kmer_index

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