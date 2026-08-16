from qc import gc_content, gc_zscores, flag_gc_outliers, gc_outlier_warnings, per_read_gc_content


def test_gc_content_all_gc():
    assert gc_content("GCGCGC") == 100.0


def test_gc_content_all_at():
    assert gc_content("ATATAT") == 0.0


def test_gc_content_mixed():
    assert gc_content("ATGC") == 50.0


def test_gc_content_handles_lowercase():
    assert gc_content("atgc") == 50.0


def test_gc_content_empty_sequence():
    assert gc_content("") == 0.0


def test_gc_zscores_known_values():
    # mean=20, sample variance=((10-20)^2+0+(30-20)^2)/(3-1)=100, std=10
    zscores = gc_zscores({"a": 10.0, "b": 20.0, "c": 30.0})
    assert zscores["a"] == -1.0
    assert zscores["b"] == 0.0
    assert zscores["c"] == 1.0


def test_gc_zscores_uniform_values_all_zero():
    zscores = gc_zscores({"a": 50.0, "b": 50.0, "c": 50.0})
    assert zscores == {"a": 0.0, "b": 0.0, "c": 0.0}


def test_gc_zscores_single_species_returns_zero():
    assert gc_zscores({"a": 50.0}) == {"a": 0.0}


def test_gc_zscores_empty_dict_returns_empty():
    assert gc_zscores({}) == {}


def test_flag_gc_outliers_flags_both_high_and_low():
    # zscores for {a:10, b:20, c:30} are -1.0, 0.0, 1.0 (see test_gc_zscores_known_values)
    flags = flag_gc_outliers({"a": 10.0, "b": 20.0, "c": 30.0}, threshold=0.5)
    assert flags == {"a": True, "b": False, "c": True}


def test_flag_gc_outliers_respects_default_threshold():
    # same zscores as above (-1.0, 0.0, 1.0) - none exceed the default 1.5 threshold
    flags = flag_gc_outliers({"a": 10.0, "b": 20.0, "c": 30.0})
    assert flags == {"a": False, "b": False, "c": False}


def test_flag_gc_outliers_uniform_values_never_flagged():
    flags = flag_gc_outliers({"a": 50.0, "b": 50.0, "c": 50.0}, threshold=0.0)
    assert flags == {"a": False, "b": False, "c": False}


def test_gc_outlier_warnings_only_includes_outliers():
    # zscores for {a:10, b:20, c:30} are -1.0, 0.0, 1.0
    warnings = gc_outlier_warnings({"a": 10.0, "b": 20.0, "c": 30.0}, threshold=0.5)
    assert set(warnings.keys()) == {"a", "c"}
    assert "b" not in warnings


def test_gc_outlier_warnings_reports_correct_direction():
    warnings = gc_outlier_warnings({"a": 10.0, "b": 20.0, "c": 30.0}, threshold=0.5)
    assert "low" in warnings["a"]
    assert "high" in warnings["c"]


def test_gc_outlier_warnings_includes_actual_numbers():
    warnings = gc_outlier_warnings({"a": 10.0, "b": 20.0, "c": 30.0}, threshold=0.5)
    assert "10.0%" in warnings["a"]
    assert "-1.00" in warnings["a"]


def test_gc_outlier_warnings_empty_when_no_outliers():
    assert gc_outlier_warnings({"a": 10.0, "b": 20.0, "c": 30.0}) == {}


def test_per_read_gc_content_maps_each_read():
    result = per_read_gc_content({"read1": "GCGC", "read2": "ATAT"})
    assert result == {"read1": 100.0, "read2": 0.0}


def test_per_read_gc_content_empty_dict_returns_empty():
    assert per_read_gc_content({}) == {}
