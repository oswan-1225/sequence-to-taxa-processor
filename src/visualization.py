import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

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

def plot_single_sample_abundance(abundance_df: pd.DataFrame, unclassified_percent: float, top_n: int = 10, title: str = "Species Abundance", output_path: str | None = None) -> Figure:
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

    Returns:
        matplotlib.figure.Figure: the chart, for the caller to display
            and/or save.
    """
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

    SPECIES_COLOR = "#2a78d6"
    UNCLASSIFIED_COLOR = "#898781"
    colors = [UNCLASSIFIED_COLOR if label == "Unclassified" else SPECIES_COLOR for label in combined.index]

    fig, ax = plt.subplots(figsize=(10, max(4, 0.4 * len(combined))))
    y_pos = range(len(combined))
    ax.barh(y_pos, list(combined.values), color=colors)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(combined.index)
    ax.invert_yaxis()
    ax.set_xlabel("Abundance (% of all reads)")
    ax.set_title(title)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    for i, value in enumerate(combined.values):
        ax.annotate(f"{value:.1f}%", (value, i), xytext=(4, 0), textcoords="offset points",
                    va="center", fontsize=9, color="#52514e")

    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, bbox_inches="tight")
    return fig
