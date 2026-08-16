import numpy as np
import pandas as pd
import pytest

from diversity import species_richness, shannon_diversity, diversity_by_sample


def test_species_richness_counts_nonzero_species():
    assert species_richness([5, 0, 3, 0]) == 2


def test_species_richness_with_threshold():
    assert species_richness([5, 2, 3], threshold=2) == 2


def test_shannon_diversity_single_species_is_zero():
    assert shannon_diversity([10]) == pytest.approx(0.0)


def test_shannon_diversity_even_split_equals_ln_n():
    assert shannon_diversity([10, 10, 10, 10]) == pytest.approx(np.log(4))


def test_shannon_diversity_empty_sample_is_zero():
    assert shannon_diversity([]) == pytest.approx(0.0)


def test_shannon_diversity_ignores_zero_count_species():
    # a species with 0 reads shouldn't change the result vs. omitting it entirely
    with_zero = shannon_diversity([10, 10, 0])
    without_zero = shannon_diversity([10, 10])
    assert with_zero == pytest.approx(without_zero)


def test_diversity_by_sample_returns_one_row_per_sample():
    abundance_df = pd.DataFrame([
        {"sample_id": "sample1", "best_match": "species_a", "count": 10, "total": 40, "percent": 25.0},
        {"sample_id": "sample1", "best_match": "species_b", "count": 10, "total": 40, "percent": 25.0},
        {"sample_id": "sample1", "best_match": "species_c", "count": 10, "total": 40, "percent": 25.0},
        {"sample_id": "sample1", "best_match": "species_d", "count": 10, "total": 40, "percent": 25.0},
        {"sample_id": "sample2", "best_match": "species_only", "count": 10, "total": 10, "percent": 100.0},
    ])

    result = diversity_by_sample(abundance_df)

    assert set(result["sample_id"]) == {"sample1", "sample2"}

    row1 = result[result["sample_id"] == "sample1"].iloc[0]
    row2 = result[result["sample_id"] == "sample2"].iloc[0]

    assert row1["species_richness"] == 4
    assert row1["shannon_diversity"] == pytest.approx(np.log(4))
    assert row2["species_richness"] == 1
    assert row2["shannon_diversity"] == pytest.approx(0.0)
