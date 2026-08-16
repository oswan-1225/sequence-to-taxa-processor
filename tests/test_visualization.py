import matplotlib
matplotlib.use("Agg")

import pandas as pd
from matplotlib.figure import Figure

from visualization import plot_species_abundance


def make_abundance_df():
    return pd.DataFrame([
        {"sample_id": "sample1", "best_match": "species_a", "count": 60, "total": 100, "percent": 60.0},
        {"sample_id": "sample1", "best_match": "species_b", "count": 40, "total": 100, "percent": 40.0},
        {"sample_id": "sample2", "best_match": "species_a", "count": 100, "total": 100, "percent": 100.0},
    ])


def test_plot_species_abundance_returns_figure():
    fig = plot_species_abundance(make_abundance_df())

    assert isinstance(fig, Figure)


def test_plot_species_abundance_saves_file_when_output_path_given(tmp_path):
    output_path = tmp_path / "abundance.png"

    plot_species_abundance(make_abundance_df(), output_path=str(output_path))

    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_plot_species_abundance_omits_other_when_top_n_covers_everything():
    fig = plot_species_abundance(make_abundance_df(), top_n=10)

    ax = fig.axes[0]
    legend_labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert "Other" not in legend_labels


def test_plot_species_abundance_collapses_excluded_species_into_other():
    df = pd.DataFrame([
        {"sample_id": "sample1", "best_match": "species_a", "count": 50, "total": 100, "percent": 50.0},
        {"sample_id": "sample1", "best_match": "species_b", "count": 30, "total": 100, "percent": 30.0},
        {"sample_id": "sample1", "best_match": "species_c", "count": 20, "total": 100, "percent": 20.0},
    ])

    fig = plot_species_abundance(df, top_n=1)

    ax = fig.axes[0]
    legend_labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert legend_labels == ["species_a", "Other"]
