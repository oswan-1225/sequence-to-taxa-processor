"""
One-off portfolio graphic: redistributed observed vs. expected per-species
abundance against the ZymoBIOMICS D6300 mock community, built from a
completed classify_reads.py run's redistributed-abundance CSV.

Dataset-specific - the 12%/2% expected values are Zymo's published
composition for this one commercial reference standard, which is why this
lives here and not in src/visualization.py (general plotting code takes
observed/expected as plain arguments; it doesn't know "Zymo" from any
other dataset). Lives in docs/ alongside the image it generates
(abundance_validation.png, embedded in README.md) rather than scripts/,
since it's meant to be re-run to refresh that specific portfolio asset,
not a disposable diagnostic.
"""
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.figure import Figure

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "src"))

REDISTRIBUTED_CSV_PATH = os.path.join(PROJECT_ROOT, "results", "ten_species_redistribution", "classified_redistributed.csv")
OUTPUT_PATH = os.path.join(PROJECT_ROOT, "results", "abundance_validation.png")

# 10 species in the ZymoBIOMICS D6300 mock community as of the 2026-08-17
# yeast reintroduction: 8 bacteria at 12% each + 2 yeasts at 2% each
# (8*12 + 2*2 = 100) - Zymo's published design, not something the
# classifier measures, a fact about the reference standard itself.
EXPECTED_ABUNDANCE = {
    "Bacillus_subtilis_complete_genome": 12.0,
    "Enterococcus_faecalis_complete_genome": 12.0,
    "Escherichia_coli_complete_genome": 12.0,
    "Lactobacillus_fermentum_complete_genome": 12.0,
    "Listeria_monocytogenes_complete_genome": 12.0,
    "Pseudomonas_aeruginosa_complete_genome": 12.0,
    "Salmonella_enterica_complete_genome": 12.0,
    "Staphylococcus_aureus_complete_genome": 12.0,
    "Cryptococcus_neoformans_complete_genome": 2.0,
    "Saccharomyces_cerevisiae_complete_genome": 2.0,
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

    # Labels sit above the dot (not left/right) so they never cross the connecting line.
    for y, o in zip(y_positions, observed_vals):
        ax.annotate(f"{o:.1f}%", (o, y), xytext=(-2, 9), textcoords="offset points",
                    va="bottom", ha="left", fontsize=9, color="#52514e")

    mean_abs_dev = sum(abs(o - e) for o, e in zip(observed_vals, expected_vals)) / len(species)
    if n_reads is not None:
        subtitle = f"{n_reads:,} Illumina sequencing reads classified — mean deviation {mean_abs_dev:.1f} percentage points"
    else:
        subtitle = f"mean deviation {mean_abs_dev:.1f} percentage points from expected"

    fig.suptitle(title, fontsize=14, fontweight="bold", y=0.98)
    ax.set_title(subtitle, fontsize=10, color="#52514e")
    ax.set_xlabel("Abundance (%)")
    ax.legend()
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, bbox_inches="tight")
    return fig


def main():
    redistributed_df = pd.read_csv(REDISTRIBUTED_CSV_PATH)
    observed = dict(zip(redistributed_df["species"], redistributed_df["percentage"]))
    # Every classified read contributes exactly 1.0 total credit across the
    # species it hit (proportional vote-share splitting), so the sum of
    # estimated_reads is exactly the total classified-read count.
    n_reads = round(redistributed_df["estimated_reads"].sum())

    plot_abundance_comparison(
        observed, EXPECTED_ABUNDANCE,
        title="Redistributed vs. Expected Abundance (ZymoBIOMICS D6300)",
        n_reads=n_reads,
        output_path=OUTPUT_PATH,
    )
    print(f"Saved to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
