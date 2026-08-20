import matplotlib
matplotlib.use("Agg")

import pandas as pd
import pytest
from matplotlib.figure import Figure

from visualization import plot_multi_sample_abundance, plot_single_sample_abundance


def make_abundance_df():
    return pd.DataFrame([
        {"sample_id": "sample1", "best_match": "species_a", "count": 60, "total": 100, "percent": 60.0},
        {"sample_id": "sample1", "best_match": "species_b", "count": 40, "total": 100, "percent": 40.0},
        {"sample_id": "sample2", "best_match": "species_a", "count": 100, "total": 100, "percent": 100.0},
    ])


def test_plot_multi_sample_abundance_returns_figure():
    fig = plot_multi_sample_abundance(make_abundance_df(), unclassified_percent={"sample1": 0.0, "sample2": 0.0})

    assert isinstance(fig, Figure)


def test_plot_multi_sample_abundance_saves_file_when_output_path_given(tmp_path):
    output_path = tmp_path / "abundance.png"

    plot_multi_sample_abundance(make_abundance_df(), unclassified_percent={"sample1": 0.0, "sample2": 0.0},
                                 output_path=str(output_path))

    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_plot_multi_sample_abundance_omits_other_when_top_n_covers_everything():
    fig = plot_multi_sample_abundance(make_abundance_df(), unclassified_percent={"sample1": 0.0, "sample2": 0.0},
                                       top_n=10)

    ax = fig.axes[0]
    legend_labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert "Other" not in legend_labels


def test_plot_multi_sample_abundance_collapses_excluded_species_into_other():
    df = pd.DataFrame([
        {"sample_id": "sample1", "best_match": "species_a", "count": 50, "total": 100, "percent": 50.0},
        {"sample_id": "sample1", "best_match": "species_b", "count": 30, "total": 100, "percent": 30.0},
        {"sample_id": "sample1", "best_match": "species_c", "count": 20, "total": 100, "percent": 20.0},
    ])

    fig = plot_multi_sample_abundance(df, unclassified_percent={"sample1": 0.0}, top_n=1)

    ax = fig.axes[0]
    legend_labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert legend_labels == ["species_a", "Other", "Unclassified"]


def test_plot_multi_sample_abundance_includes_unclassified_segment():
    df = pd.DataFrame([
        {"sample_id": "sample1", "best_match": "species_a", "count": 80, "total": 80, "percent": 100.0},
    ])

    fig = plot_multi_sample_abundance(df, unclassified_percent={"sample1": 20.0})

    ax = fig.axes[0]
    legend_labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert "Unclassified" in legend_labels

    bar_container = ax.containers[-1]
    assert bar_container.get_label() == "Unclassified"
    assert bar_container[0].get_height() == pytest.approx(20.0)


def make_single_sample_df():
    return pd.DataFrame([
        {"sample_id": "sample1", "best_match": "species_a", "count": 60, "total": 80, "percent": 75.0},
        {"sample_id": "sample1", "best_match": "species_b", "count": 20, "total": 80, "percent": 25.0},
    ])


def bar_widths_by_label(fig):
    ax = fig.axes[0]
    labels = [t.get_text() for t in ax.get_yticklabels()]
    widths = [p.get_width() for p in ax.patches]
    return dict(zip(labels, widths))


def test_plot_single_sample_abundance_returns_figure():
    fig = plot_single_sample_abundance(make_single_sample_df(), unclassified_percent=10.0)

    assert isinstance(fig, Figure)


def test_plot_single_sample_abundance_saves_file_when_output_path_given(tmp_path):
    output_path = tmp_path / "abundance.png"

    plot_single_sample_abundance(make_single_sample_df(), unclassified_percent=10.0, output_path=str(output_path))

    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_plot_single_sample_abundance_includes_unclassified_bar():
    fig = plot_single_sample_abundance(make_single_sample_df(), unclassified_percent=10.0)

    labels = [t.get_text() for t in fig.axes[0].get_yticklabels()]
    assert "Unclassified" in labels


def test_plot_single_sample_abundance_rescales_to_percent_of_all_reads():
    fig = plot_single_sample_abundance(make_single_sample_df(), unclassified_percent=10.0)

    widths = bar_widths_by_label(fig)
    assert widths["species_a"] == pytest.approx(67.5)  # 75% of classified * 0.9 classified fraction
    assert widths["species_b"] == pytest.approx(22.5)  # 25% of classified * 0.9 classified fraction
    assert widths["Unclassified"] == pytest.approx(10.0)


def test_plot_single_sample_abundance_collapses_excluded_species_into_other():
    fig = plot_single_sample_abundance(make_single_sample_df(), unclassified_percent=10.0, top_n=1)

    labels = [t.get_text() for t in fig.axes[0].get_yticklabels()]
    assert "Other" in labels
    assert "species_b" not in labels


def test_plot_single_sample_abundance_omits_other_when_top_n_covers_everything():
    fig = plot_single_sample_abundance(make_single_sample_df(), unclassified_percent=10.0, top_n=10)

    labels = [t.get_text() for t in fig.axes[0].get_yticklabels()]
    assert "Other" not in labels
