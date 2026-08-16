import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

def plot_species_abundance(abundance_df: pd.DataFrame, top_n: int = 10, title: str = "Species Abundance", output_path: str | None = None) -> Figure:
    """
    Stacked bar chart of per-species relative abundance (%) for each sample.

    Species outside the top_n (by total abundance across samples) are
    collapsed into a single "Other" row, unless every excluded species is
    at 0% (in which case "Other" is omitted rather than plotted as an
    empty segment).

    Parameters:
        abundance_df: DataFrame from database.get_abundance(), with columns
            ['sample_id', 'best_match', 'count', 'total', 'percent'] - one
            row per (sample, species).
        top_n: number of top species to display individually; the rest are
            summed into "Other".
        title: title of the plot.
        output_path: if given, save the figure to this path (e.g. .png)
            before returning it.

    Returns:
        matplotlib.figure.Figure: the chart, for the caller to display
            and/or save.
    """

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

    fig, ax = plt.subplots(figsize=(12, 8))
    top_pivot.T.plot(kind='bar', stacked=True, ax=ax)
    ax.set_title(title)
    ax.set_xlabel('Sample')
    ax.set_ylabel('Abundance (%)')
    ax.legend(title='Species', bbox_to_anchor=(1.02, 1), loc='upper left')
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, bbox_inches='tight')
    return fig
