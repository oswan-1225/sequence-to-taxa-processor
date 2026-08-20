import math
import textwrap
from typing import cast

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

SPECIES_COLOR = "#2a78d6"
UNCLASSIFIED_COLOR = "#898781"
INK = "#2b2a28"
MUTED = "#52514e"


def _gc_flagged_species(abundance_df: pd.DataFrame) -> dict[str, str]:
    """
    Species with a non-empty GC-outlier warning, mapped to the warning text.

    qc.gc_outlier_warnings() deliberately omits non-outliers, so pipeline.py's
    .map() leaves NaN in the 'gc_warning' column for every species that passed.
    Returns {} when the column is absent entirely (a caller that never ran the
    GC check).

    Parameters:
        abundance_df: DataFrame from database.get_abundance(), optionally with
            a 'gc_warning' column added by run_pipeline().

    Returns:
        dict[str, str]: species name -> warning text, for flagged species only.
    """
    if "gc_warning" not in abundance_df.columns:
        return {}

    flagged = abundance_df.dropna(subset=["gc_warning"])
    flagged = flagged[flagged["gc_warning"].astype(str).str.strip() != ""]
    return dict(zip(flagged["best_match"], flagged["gc_warning"].astype(str)))


def _gc_footnote_text(flagged: dict[str, str], width: int = 118) -> str | None:
    """
    Wrapped footnote text for the GC-outlier warnings marked on the bars.

    Parameters:
        flagged: species -> warning text, from _gc_flagged_species().
        width: characters per line before wrapping.

    Returns:
        str | None: the footnote, or None if nothing was flagged.
    """
    if not flagged:
        return None

    lines = []
    for species, message in flagged.items():
        lines.append(textwrap.fill(f"* {species}: {message}", width=width,
                                    subsequent_indent="  "))
    return "\n".join(lines)


def _footnote_height_inches(footnote: str | None) -> float:
    """
    Vertical space a footnote needs, in inches, for figure sizing.

    Parameters:
        footnote: text from _gc_footnote_text(), or None.

    Returns:
        float: inches to add to the figure height; 0.0 when there is no
            footnote.
    """
    if not footnote:
        return 0.0
    return 0.2 * len(footnote.splitlines()) + 0.15


def _layout_with_footnote(fig: Figure, footnote: str | None, reserved: float) -> None:
    """
    Lay out a figure, keeping the bottom clear for a footnote and drawing it.

    Callers size the figure to include _footnote_height_inches() already, so
    this only has to stop tight_layout from reclaiming that strip.

    Parameters:
        fig: the figure to lay out.
        footnote: text to draw at the bottom left, or None for a plain layout.
        reserved: fraction of the figure height to keep clear at the bottom.
    """
    if not footnote:
        fig.tight_layout()
        return

    fig.text(0.01, 0.01, footnote, fontsize=7.5, color=MUTED, va="bottom", ha="left")
    fig.tight_layout(rect=(0, reserved, 1, 1))


def plot_multi_sample_abundance(abundance_df: pd.DataFrame, unclassified_percent: dict[str, float] | pd.Series,
                                 top_n: int = 10, title: str = "Species Abundance", output_path: str | None = None) -> Figure:
    """
    Stacked bar chart of per-species relative abundance (% of all reads,
    including unclassified) for each sample. The multi-sample sibling of
    plot_single_sample_abundance - use that one instead for a single
    sample, which pipeline.py always is today.

    Species outside the top_n (by total abundance across samples) are
    collapsed into a single "Other" row, unless every excluded species is
    at 0% (in which case "Other" is omitted rather than plotted as an
    empty segment). Every sample also gets an "Unclassified" segment.

    Parameters:
        abundance_df: DataFrame from database.get_abundance(), with columns
            ['sample_id', 'best_match', 'count', 'total', 'percent'] - one
            row per (sample, species). 'percent' is classified-reads-only,
            rescaled here against unclassified_percent the same way
            plot_single_sample_abundance does.
        unclassified_percent: sample_id -> unclassified share (% of all
            reads), from database.get_classification_totals().
        top_n: number of top species to display individually; the rest are
            summed into "Other".
        title: title of the plot.
        output_path: if given, save the figure to this path (e.g. .png)
            before returning it.

    Returns:
        matplotlib.figure.Figure: the chart, for the caller to display
            and/or save.
    """
    abundance_df = abundance_df.copy()
    scale = 1 - (abundance_df['sample_id'].map(unclassified_percent) / 100)
    abundance_df['percent'] = abundance_df['percent'] * scale

    pivot = abundance_df.pivot_table(index='best_match', columns='sample_id', values='percent')
    pivot.fillna(0.0, inplace=True)

    totals = pivot.sum(axis=1)
    top_species = totals.nlargest(top_n).index
    top_pivot = pivot.loc[top_species]
    excluded = pivot.copy()
    excluded.drop(index=top_species, inplace=True)
    total_others = excluded.sum(axis=0)

    if total_others.any():
        top_pivot.loc['Other'] = total_others

    top_pivot.loc['Unclassified'] = pd.Series(unclassified_percent).reindex(pivot.columns)

    fig, ax = plt.subplots(figsize=(12, 8))
    top_pivot.T.plot(kind='bar', stacked=True, ax=ax)
    ax.set_title(title)
    ax.set_xlabel('Sample')
    ax.set_ylabel('Abundance (% of all reads)')
    ax.legend(title='Species', bbox_to_anchor=(1.02, 1), loc='upper left')
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, bbox_inches='tight')
    return fig

def plot_single_sample_abundance(abundance_df: pd.DataFrame, unclassified_percent: float, top_n: int = 10, title: str = "Species Abundance", output_path: str | None = None, ax: Axes | None = None) -> Figure:
    """
    Horizontal bar chart of per-species abundance for a single sample,
    ranked highest to lowest, with an explicit "Unclassified" bar.

    Every bar (species, "Other", "Unclassified") is expressed as a share
    of the sample's total reads, abundance_df's own 'percent' column is classified-reads-only, so it's
    rescaled here to be comparable with unclassified_percent (already a
    share of all reads).

    Parameters:
        abundance_df: DataFrame from database.get_abundance() for a single
            sample - columns ['sample_id', 'best_match', 'count', 'total',
            'percent'].
        unclassified_percent: this sample's unclassified share (% of all
            reads), from database.get_classification_totals().
        top_n: number of species to show individually; the rest collapse
            into "Other".
        title: title of the plot.
        output_path: if given, save the figure to this path (e.g. .png)
            before returning it.
        ax: draw into this existing axes instead of creating a new figure.
            When given, the caller owns the figure - layout and any GC
            footnote are left to it (see plot_sample_summary, which composes
            this chart under a row of stat tiles). When None, this function
            makes its own figure and adds the footnote itself.

    If abundance_df carries a 'gc_warning' column (run_pipeline() adds one),
    flagged species get a "*" appended to their bar label, keyed to a footnote
    naming the species and the reason.

    Returns:
        matplotlib.figure.Figure: the chart, for the caller to display
            and/or save.
    """
    flagged = _gc_flagged_species(abundance_df)

    abundance_df = abundance_df.copy()
    scale = 1 - (unclassified_percent / 100)
    abundance_df["percent_of_total_reads"] = abundance_df["percent"] * scale

    species_pct = abundance_df.set_index("best_match")["percent_of_total_reads"]
    top_species = species_pct.nlargest(top_n)
    other = species_pct.drop(top_species.index)

    combined = top_species.copy()
    if other.sum() > 0:
        combined["Other"] = other.sum()
    combined["Unclassified"] = unclassified_percent
    combined = combined.sort_values(ascending=False)

    colors = [UNCLASSIFIED_COLOR if label == "Unclassified" else SPECIES_COLOR for label in combined.index]

    owns_figure = ax is None
    footnote = _gc_footnote_text(flagged) if owns_figure else None
    bar_h = max(4.0, 0.4 * len(combined))
    foot_h = _footnote_height_inches(footnote)
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, bar_h + foot_h))
    else:
        fig = cast(Figure, ax.figure)

    y_pos = range(len(combined))
    ax.barh(y_pos, list(combined.values), color=colors)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(combined.index)
    ax.invert_yaxis()
    ax.set_xlabel("Abundance (% of all reads)")
    ax.set_title(title)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    for i, (label, value) in enumerate(zip(combined.index, combined.values)):
        marker = " *" if label in flagged else ""
        ax.annotate(f"{value:.1f}%{marker}", (value, i), xytext=(4, 0), textcoords="offset points",
                    va="center", fontsize=9, color=MUTED)

    if owns_figure:
        _layout_with_footnote(fig, footnote, foot_h / (bar_h + foot_h))

    if output_path:
        fig.savefig(output_path, bbox_inches="tight")
    return fig


def _draw_stat_tiles(ax: Axes, tiles: list[tuple[str, str, str | None]],
                     pitch: float = 0.25) -> None:
    """
    Render a row of stat tiles across a blank axes.

    Each tile is (label, value, sub): a muted sentence-case label, a large
    semibold value, and an optional smaller line underneath. Text uses the
    ink/muted grays only - never a series color - so the tiles never read as
    another data encoding.

    Parameters:
        ax: a dedicated axes; its frame and ticks are turned off here.
        tiles: the tiles to draw, left to right.
        pitch: horizontal spacing between tile origins, as a fraction of the
            axes width. Fixed rather than 1/len(tiles) so a row of two or
            three tiles stays grouped at the left instead of being stretched
            across the full width.
    """
    ax.set_axis_off()
    if not tiles:
        return

    for i, (label, value, sub) in enumerate(tiles):
        x = i * pitch
        ax.text(x, 0.92, label, transform=ax.transAxes, fontsize=9.5,
                color=MUTED, va="top", ha="left")
        ax.text(x, 0.66, value, transform=ax.transAxes, fontsize=21,
                color=INK, va="top", ha="left", fontweight="semibold")
        if sub:
            ax.text(x, 0.04, sub, transform=ax.transAxes, fontsize=8.5,
                    color=MUTED, va="bottom", ha="left")


def plot_sample_summary(abundance_df: pd.DataFrame, stats: pd.Series | dict, top_n: int = 10,
                        title: str = "Species Abundance", output_path: str | None = None) -> Figure:
    """
    Single-sample abundance chart with a row of summary stat tiles above it.

    Composes plot_single_sample_abundance() under three tiles carrying the
    numbers the bars themselves cannot show:

      - Reads classified: the classification rate plus raw sample size. Every
        bar is a percentage, so nothing else on the figure says how many reads
        the percentages rest on.
      - Species detected: species_richness. Distinct from the bar count
        whenever richness exceeds top_n, since the excess collapses into
        "Other" - the sub-line says so when that happens.
      - Evenness: Pielou's J = H / ln(S), the share of the maximum Shannon
        diversity achievable with S species. Raw H conflates richness with
        evenness and is unitless; dividing by ln(S) removes the richness
        term and bounds the result to 0-1, so it compares across samples with
        different species counts. H is kept as the sub-line since it is what
        diversity_report.csv reports.

    Unclassified percent deliberately gets no tile - the chart already has an
    explicit "Unclassified" bar, and the rate is folded into the first tile.

    Parameters:
        abundance_df: DataFrame from database.get_abundance() for a single
            sample, optionally with run_pipeline()'s added 'gc_warning'
            column (which drives the bar markers and footnote).
        stats: one sample's row from the merged diversity report - a Series or
            dict with keys 'total_reads', 'classified_reads',
            'unclassified_percent', 'species_richness', 'shannon_diversity'.
        top_n: number of species to plot individually before collapsing the
            rest into "Other".
        title: title of the abundance chart beneath the tiles.
        output_path: if given, save the figure to this path before returning it.

    Returns:
        matplotlib.figure.Figure: the composed figure.
    """
    total_reads = int(stats["total_reads"])
    classified_reads = int(stats["classified_reads"])
    unclassified_percent = float(stats["unclassified_percent"])
    richness = int(stats["species_richness"])
    shannon = float(stats["shannon_diversity"])

    if total_reads > 0:
        rate = f"{classified_reads / total_reads * 100:.1f}%"
    else:
        rate = "n/a"

    # ln(1) is 0, so evenness is undefined for a single species (and for none).
    evenness = f"{shannon / math.log(richness):.2f}" if richness > 1 else "n/a"

    tiles: list[tuple[str, str, str | None]] = [
        ("Reads classified", rate, f"{classified_reads:,} of {total_reads:,}"),
        ("Species detected", f"{richness:,}",
         f"showing top {top_n} individually" if richness > top_n else None),
        ("Evenness (Pielou J)", evenness, f"Shannon H = {shannon:.2f}"),
    ]

    n_species_shown = min(len(abundance_df), top_n)
    n_bars = n_species_shown + (1 if len(abundance_df) > top_n else 0) + 1

    flagged = _gc_flagged_species(abundance_df)
    footnote = _gc_footnote_text(flagged)

    tile_h = 1.15
    bar_h = max(4.0, 0.4 * n_bars)
    foot_h = _footnote_height_inches(footnote)

    fig, (tile_ax, bar_ax) = plt.subplots(
        2, 1,
        figsize=(10, tile_h + bar_h + foot_h),
        gridspec_kw={"height_ratios": [tile_h, bar_h]},
    )

    _draw_stat_tiles(tile_ax, tiles)
    plot_single_sample_abundance(abundance_df, unclassified_percent, top_n=top_n,
                                  title=title, ax=bar_ax)

    _layout_with_footnote(fig, footnote, foot_h / (tile_h + bar_h + foot_h))

    # tight_layout aligns the tile axes with the bar axes, which is inset by
    # the width of the species labels. Widen it back out afterwards so the
    # tiles start at the figure's left edge instead of floating mid-figure.
    tile_pos = tile_ax.get_position()
    tile_ax.set_position((0.015, tile_pos.y0, 0.97, tile_pos.height))

    if output_path:
        fig.savefig(output_path, bbox_inches="tight")
    return fig
