"""
One-off portfolio graphic: winner-take-all observed vs. expected per-species
abundance against the ZymoBIOMICS D6300 mock community, built from a
completed classify_reads.py run's classifications CSV.

Dataset-specific - the expected values are Zymo's own published "16S Only"
composition (rRNA copy number-adjusted per species, not a flat share),
which is why this lives here and not in src/visualization.py (general
plotting code takes observed/expected as plain arguments; it doesn't know
"Zymo" from any other dataset). SRR10391187 is confirmed 16S rRNA amplicon
sequencing, so this is the theoretical composition Zymo's own
documentation says to use, not the flat genomic-DNA composition an earlier
version of this script used. The 2 yeasts are excluded: Zymo's 16S-only
column lists them as NA (bacteria-targeted 16S primers don't amplify
fungal rRNA), so there's no valid expected value to plot them against.
Source: https://files.zymoresearch.com/protocols/_d6300_zymobiomics_microbial_community_standard.pdf
(Table 1, "16S Only" column). Lives in docs/ alongside the image it
generates (abundance_validation.png, embedded in README.md) rather than
scripts/, since it's meant to be re-run to refresh that specific portfolio
asset, not a disposable diagnostic.
"""
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.figure import Figure

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "src"))

# Default run folder this graphic is built from. Override on the command line
# to plot a different run: python docs/plot_zymobiomics_validation.py <csv>
DEFAULT_CLASSIFIED_CSV = os.path.join(PROJECT_ROOT, "results", "zymo_run_a", "classifications.csv")
OUTPUT_PATH = os.path.join(PROJECT_ROOT, "results", "abundance_validation.png")

# Zymo's published "16S Only" theoretical composition for D6300 (Table 1),
# rRNA operon copy number-adjusted per species - the correct expected
# values for 16S targeted sequencing, per Zymo's own footnote on that
# column. Sums to 100% across these 8 bacteria alone.
EXPECTED_ABUNDANCE = {
    "Pseudomonas_aeruginosa_complete_genome": 4.2,
    "Escherichia_coli_complete_genome": 10.1,
    "Salmonella_enterica_complete_genome": 10.4,
    "Lactobacillus_fermentum_complete_genome": 18.4,
    "Enterococcus_faecalis_complete_genome": 9.9,
    "Staphylococcus_aureus_complete_genome": 15.5,
    "Listeria_monocytogenes_complete_genome": 14.1,
    "Bacillus_subtilis_complete_genome": 17.4,
}


def _display_name(species: str) -> str:
    """
    'Pseudomonas_aeruginosa_complete_genome' -> 'Pseudomonas aeruginosa'.
    Strips the '_complete_genome' suffix these reference genome filenames
    all share and swaps underscores for spaces - a fact about this
    dataset's file naming, not something src/visualization.py should know.
    """
    return species.replace("_complete_genome", "").replace("_", " ")


def plot_abundance_comparison(observed: dict[str, float], expected: dict[str, float],
                               title: str = "Observed vs Expected Abundance",
                               n_reads: int | None = None,
                               output_path: str | None = None) -> Figure:
    """
    Dumbbell chart comparing each species' observed classified abundance (%)
    against its expected abundance (%), sorted by the size of the miss
    (largest deviation first, so the reader's eye lands there).
    """
    species = sorted(observed.keys(), key=lambda s: abs(observed[s] - expected[s]), reverse=True)
    expected_vals = [expected[s] for s in species]
    observed_vals = [observed[s] for s in species]

    fig, ax = plt.subplots(figsize=(10, 6))
    y_positions = list(range(len(species)))

    ax.set_yticks(y_positions)
    ax.set_yticklabels([_display_name(s) for s in species])
    ax.invert_yaxis()
    # Extra top/bottom margin so the first/last row's labels have room
    # above/below the dot instead of clipping against the x-axis.
    ax.set_ylim(y_positions[-1] + 0.7, y_positions[0] - 0.7)

    # Zebra striping and gridlines behind everything else.
    for y in y_positions[::2]:
        ax.axhspan(y - 0.5, y + 0.5, color="#c3c2b7", alpha=0.12, zorder=0)
    ax.xaxis.grid(True, color="#e8e7e0", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)

    xmin = [min(e, o) for e, o in zip(expected_vals, observed_vals)]
    xmax = [max(e, o) for e, o in zip(expected_vals, observed_vals)]

    ax.hlines(y_positions, xmin, xmax, colors="#c3c2b7", linewidth=2, zorder=1)
    # White ring keeps each dot legible crossing the connecting line/gridlines.
    ax.scatter(expected_vals, y_positions, color="#4B4DA0", s=100,
               edgecolors="white", linewidths=1.5, zorder=3, label="Expected")
    ax.scatter(observed_vals, y_positions, color="#5F9930", s=60,
               edgecolors="white", linewidths=1.5, zorder=3, label="Observed")

    # Start the x-axis at 0 rather than autoscaling tight to the data - a
    # truncated axis makes any gap look more dramatic than it is.
    xmax_val = max(expected_vals + observed_vals) * 1.2
    ax.set_xlim(0, xmax_val)

    # Both dots get direct labels now that Expected varies per species
    # instead of repeating one flat number. Observed is the actual
    # measurement, so it's bold and full-size, centered above its dot.
    # Expected is reference/context, so it's italic, smaller, and lighter,
    # centered below its own dot - a real size/weight hierarchy instead of
    # two equally-loud numbers competing for attention.
    for y, o in zip(y_positions, observed_vals):
        ax.annotate(f"{o:.1f}%", (o, y), xytext=(0, 9), textcoords="offset points",
                    va="bottom", ha="center", fontsize=9, fontweight="bold", color="#52514e")
    for y, e in zip(y_positions, expected_vals):
        ax.annotate(f"{e:.1f}%", (e, y), xytext=(0, -8), textcoords="offset points",
                    va="top", ha="center", fontsize=8, fontstyle="italic", color="#8a8975")

    mean_rel_dev = sum(abs(o - e) / e for o, e in zip(observed_vals, expected_vals)) / len(species) * 100
    if n_reads is not None:
        subtitle = f"{n_reads:,} Illumina sequencing reads classified, mean relative deviation {mean_rel_dev:.1f}%"
    else:
        subtitle = f"mean relative deviation {mean_rel_dev:.1f}% from expected"

    fig.suptitle(title, fontsize=14, fontweight="bold", y=0.98)
    ax.set_title(subtitle, fontsize=10, color="#52514e")
    ax.set_xlabel("Abundance (%)")
    legend = ax.legend()
    # Match the legend text styling to the direct labels (bold Observed,
    # italic Expected) so the two reinforce each other instead of the
    # legend using a third, unstyled convention of its own.
    expected_text, observed_text = legend.get_texts()
    expected_text.set_fontstyle("italic")
    observed_text.set_fontweight("bold")
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, bbox_inches="tight")
    return fig


def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_CLASSIFIED_CSV
    if not os.path.exists(csv_path):
        raise SystemExit(
            f"No classification CSV at {csv_path}. Run the pipeline first, or pass "
            f"a path: python docs/plot_zymobiomics_validation.py <classifications.csv>"
        )
    classified_df = pd.read_csv(csv_path)
    classified_df = classified_df[classified_df["best_match"].notna()]
    counts = classified_df["best_match"].value_counts()
    n_reads = int(counts.sum())
    observed = {species: count / n_reads * 100 for species, count in counts.items()
                if species in EXPECTED_ABUNDANCE}

    plot_abundance_comparison(
        observed, EXPECTED_ABUNDANCE,
        title="Observed vs. Zymo's 16S-Adjusted Expected Abundance (ZymoBIOMICS D6300)",
        n_reads=n_reads,
        output_path=OUTPUT_PATH,
    )
    print(f"Saved to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
