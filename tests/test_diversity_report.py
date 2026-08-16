import numpy as np
import pytest

from database import create_database, insert_sample_results
from diversity_report import generate_report


def test_generate_report_matches_hand_computed_diversity(tmp_path):
    db_path = str(tmp_path / "test.db")
    create_database(db_path)

    sample1_results = [
        {"read_id": "r1", "best_match": "species_a", "confidence": 1.0},
        {"read_id": "r2", "best_match": "species_a", "confidence": 1.0},
        {"read_id": "r3", "best_match": "species_b", "confidence": 1.0},
        {"read_id": "r4", "best_match": None, "confidence": 0.0},
    ]
    sample2_results = [
        {"read_id": "r5", "best_match": "species_only", "confidence": 1.0},
    ]
    insert_sample_results(db_path, "sample1", "test", sample1_results)
    insert_sample_results(db_path, "sample2", "test", sample2_results)

    report_df = generate_report(db_path)

    assert set(report_df["sample_id"]) == {"sample1", "sample2"}

    row1 = report_df[report_df["sample_id"] == "sample1"].iloc[0]
    row2 = report_df[report_df["sample_id"] == "sample2"].iloc[0]

    # sample1: 2 species_a, 1 species_b (unclassified read excluded)
    assert row1["species_richness"] == 2
    p = np.array([2 / 3, 1 / 3])
    assert row1["shannon_diversity"] == pytest.approx(-np.sum(p * np.log(p)))

    # sample2: single species -> zero diversity
    assert row2["species_richness"] == 1
    assert row2["shannon_diversity"] == pytest.approx(0.0)
