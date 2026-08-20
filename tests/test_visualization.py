import math

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import pytest
from matplotlib.figure import Figure

from visualization import (_gc_flagged_species, _gc_footnote_text, plot_multi_sample_abundance,
                           plot_sample_summary, plot_single_sample_abundance)


@pytest.fixture(autouse=True)
def close_figures():
    """pyplot keeps every figure alive until closed; without this the suite
    trips matplotlib's 20-open-figure warning partway through."""
    yield
    plt.close("all")


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


# Stand-in for whatever qc.gc_outlier_warnings() produces. These tests assert on
# the marker/footnote plumbing, never on the wording or any measured value, so
# nothing here is tied to a particular reference set.
GC_WARNING = "GC content is unusually high for this reference set."


def make_summary_stats(**overrides):
    """One sample's row of the merged diversity report, as pipeline.py passes it."""
    stats = {
        "sample_id": "sample1",
        "total_reads": 100,
        "classified_reads": 80,
        "unclassified_reads": 20,
        "unclassified_percent": 20.0,
        "species_richness": 2,
        "shannon_diversity": math.log(2),  # perfectly even across 2 species -> J = 1.00
    }
    stats.update(overrides)
    return pd.Series(stats)


def tile_texts(fig):
    """Every string drawn on the stat-tile axes (always the first axes)."""
    return [t.get_text() for t in fig.axes[0].texts]


def bar_annotations(fig):
    """Every string drawn on the composed figure's bar axes (always the second)."""
    return [t.get_text() for t in fig.axes[1].texts]


def test_plot_sample_summary_returns_figure():
    fig = plot_sample_summary(make_single_sample_df(), make_summary_stats())

    assert isinstance(fig, Figure)


def test_plot_sample_summary_saves_file_when_output_path_given(tmp_path):
    output_path = tmp_path / "summary.png"

    plot_sample_summary(make_single_sample_df(), make_summary_stats(), output_path=str(output_path))

    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_plot_sample_summary_tiles_report_classification_rate_and_sample_size():
    fig = plot_sample_summary(make_single_sample_df(), make_summary_stats())

    texts = tile_texts(fig)
    assert "Reads classified" in texts
    assert "80.0%" in texts
    assert "80 of 100" in texts


def test_plot_sample_summary_tiles_report_evenness_and_shannon():
    fig = plot_sample_summary(make_single_sample_df(), make_summary_stats())

    texts = tile_texts(fig)
    assert "1.00" in texts  # J = ln(2) / ln(2), a perfectly even 2-species sample
    assert "Shannon H = 0.69" in texts


def test_plot_sample_summary_evenness_is_not_available_for_a_single_species():
    df = pd.DataFrame([
        {"sample_id": "sample1", "best_match": "species_a", "count": 80, "total": 80, "percent": 100.0},
    ])

    fig = plot_sample_summary(df, make_summary_stats(species_richness=1, shannon_diversity=0.0))

    # ln(1) is 0, so J would divide by zero rather than mean anything.
    assert "n/a" in tile_texts(fig)


def test_plot_sample_summary_classification_rate_is_not_available_without_reads():
    df = pd.DataFrame(
        [], columns=["sample_id", "best_match", "count", "total", "percent"]
    ).astype({"count": "int64", "total": "int64", "percent": "float64"})

    fig = plot_sample_summary(df, make_summary_stats(
        total_reads=0, classified_reads=0, unclassified_reads=0,
        unclassified_percent=0.0, species_richness=0, shannon_diversity=0.0))

    assert "n/a" in tile_texts(fig)


def test_plot_sample_summary_notes_when_richness_exceeds_the_plotted_species():
    fig = plot_sample_summary(make_single_sample_df(), make_summary_stats(species_richness=50), top_n=10)

    assert "showing top 10 individually" in tile_texts(fig)


def test_plot_sample_summary_omits_the_top_n_note_when_every_species_is_plotted():
    fig = plot_sample_summary(make_single_sample_df(), make_summary_stats(species_richness=2), top_n=10)

    assert "showing top 10 individually" not in tile_texts(fig)


def test_plot_sample_summary_marks_gc_flagged_species_on_its_bar():
    df = make_single_sample_df()
    df["gc_warning"] = [GC_WARNING, None]

    fig = plot_sample_summary(df, make_summary_stats())

    marked = [t for t in bar_annotations(fig) if t.endswith(" *")]
    assert len(marked) == 1
    assert marked[0].startswith("60.0%")  # species_a, the flagged one


def test_plot_sample_summary_leaves_bars_unmarked_when_no_species_is_flagged():
    df = make_single_sample_df()
    df["gc_warning"] = [None, None]

    fig = plot_sample_summary(df, make_summary_stats())

    assert not any(t.endswith(" *") for t in bar_annotations(fig))


def test_plot_sample_summary_works_without_a_gc_warning_column():
    fig = plot_sample_summary(make_single_sample_df(), make_summary_stats())

    assert isinstance(fig, Figure)
    assert not any(t.endswith(" *") for t in bar_annotations(fig))


def test_plot_single_sample_abundance_draws_into_a_supplied_axes():
    outer, ax = plt.subplots()
    fig = plot_single_sample_abundance(make_single_sample_df(), unclassified_percent=10.0, ax=ax)

    assert fig is outer
    assert len(outer.axes) == 1


def test_gc_flagged_species_ignores_unflagged_and_missing_columns():
    df = make_single_sample_df()
    assert _gc_flagged_species(df) == {}

    df["gc_warning"] = [GC_WARNING, None]
    assert _gc_flagged_species(df) == {"species_a": GC_WARNING}


def test_gc_footnote_text_is_none_when_nothing_is_flagged():
    assert _gc_footnote_text({}) is None
    assert "species_a" in _gc_footnote_text({"species_a": GC_WARNING})
