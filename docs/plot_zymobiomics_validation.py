"""
One-off portfolio graphic: observed vs. expected per-species abundance
against the ZymoBIOMICS D6300 mock community, built from a completed
pipeline.py run's classifications database.

Dataset-specific - the ~12.5%-each expected values are Zymo's published
composition for this one commercial reference standard, which is why this
lives here and not in src/visualization.py (general plotting code takes
observed/expected as plain arguments; it doesn't know "Zymo" from any
other dataset). See scripts/investigate_pseudomonas_gc.py for the same
pattern.
"""
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.figure import Figure

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "src"))

from database import get_abundance  # type: ignore

DB_PATH = os.path.join(PROJECT_ROOT, "results", "pseudomonas_investigation", "classifications.db")
OUTPUT_PATH = os.path.join(PROJECT_ROOT, "results", "abundance_validation.png")

# 8 species in the ZymoBIOMICS D6300 mock community, evenly represented by
# Zymo's published design (1/8 = 12.5% each) - not something the classifier
# measures, a fact about the reference standard itself.
EXPECTED_ABUNDANCE = {
    "Bacillus_subtilis_complete_genome": 12.5,
    "Enterococcus_faecalis_complete_genome": 12.5,
    "Escherichia_coli_complete_genome": 12.5,
    "Lactobacillus_fermentum_complete_genome": 12.5,
    "Listeria_monocytogenes_complete_genome": 12.5,
    "Pseudomonas_aeruginosa_complete_genome": 12.5,
    "Salmonella_enterica_complete_genome": 12.5,
    "Staphylococcus_aureus_complete_genome": 12.5,
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

    xmin = [min(e, o) for e, o in zip(expected_vals, observed_vals)]
    xmax = [max(e, o) for e, o in zip(expected_vals, observed_vals)]

    ax.hlines(y_positions, xmin, xmax, colors="#c3c2b7", linewidth=2, zorder=1)
    ax.scatter(expected_vals, y_positions, color="#4B4DA0", s=100, zorder=2, label="Expected")
    ax.scatter(observed_vals, y_positions, color="#5F9930", s=100, zorder=2, label="Observed")

    # Direct labels on the Observed dots only - Expected repeats the same
    # 12.5% eight times, so labeling it too would just be noise. Offset in
    # points (screen space), not data units, so the gap next to the dot
    # stays a constant size regardless of the x-axis scale. Which side the
    # label sits on flips with which side the dot is on, so it never
    # collides with the connecting line.
    for y, e, o in zip(y_positions, expected_vals, observed_vals):
        offset, ha = (8, "left") if o >= e else (-8, "right")
        ax.annotate(f"{o:.1f}%", (o, y), xytext=(offset, 0), textcoords="offset points",
                    va="center", ha=ha, fontsize=9, color="#52514e")

    # Start the x-axis at 0 rather than autoscaling tight to the data - a
    # truncated axis makes any gap look more dramatic than it is. Extra
    # headroom (not just 10%) leaves room for the direct labels above.
    ax.set_xlim(0, max(expected_vals + observed_vals) * 1.2)

    mean_abs_dev = sum(abs(o - e) for o, e in zip(observed_vals, expected_vals)) / len(species)
    if n_reads is not None:
        subtitle = f"{n_reads:,} real Illumina sequencing reads classified — mean deviation {mean_abs_dev:.1f} percentage points"
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
    abundance_df = get_abundance(DB_PATH)
    observed = dict(zip(abundance_df["best_match"], abundance_df["percent"]))
    n_reads = int(abundance_df["total"].iloc[0])

    plot_abundance_comparison(
        observed, EXPECTED_ABUNDANCE,
        title="Classified vs. Expected Abundance (ZymoBIOMICS D6300)",
        n_reads=n_reads,
        output_path=OUTPUT_PATH,
    )
    print(f"Saved to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
