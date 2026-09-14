"""
Bootstrap confidence intervals for the three read-classification accuracy
methods measured against the ZymoBIOMICS D6300 mock community (winner-take-all,
proportional redistribution, discard-ambiguous). Turns the point estimates
already reported in the README into intervals, so a change in the number
(e.g. 14.8% -> 17.3% after canonicalizing k-mers) can be judged against
sampling noise rather than compared as if it were exact.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "src"))

DEFAULT_CLASSIFIED_CSV = os.path.join(PROJECT_ROOT, "results", "zymo_run_canonical_1", "classifications.csv")

# Zymo's published "16S Only" theoretical composition for D6300 (Table 1),
# rRNA operon copy number-adjusted per species - same values and rationale
# as plot_zymobiomics_validation.py.
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

N_BOOTSTRAP = 10_000
SEED = 42  # arbitrary but fixed, so the interval is reproducible run to run


def mean_relative_deviation(counts: pd.Series, expected: dict) -> float:
    """
    counts: species -> read count. Must include every species that could
    appear (not just the 8 scored ones) - see the note in main() about why.

    Returns the same statistic plot_zymobiomics_validation.py's subtitle
    reports: mean(|observed% - expected%| / expected%) over the 8 species
    with a defined expected value, as a percentage.
    """
    n = counts.sum()
    deviations = []
    for species, expected_pct in expected.items():
        observed_pct = counts.get(species, 0) / n * 100
        deviations.append(abs(observed_pct - expected_pct) / expected_pct)
    return sum(deviations) / len(deviations) * 100


def bootstrap_ci(counts: pd.Series, expected: dict, n_bootstrap: int,
                  rng: np.random.Generator) -> tuple[float, float, float]:
    """
    Parametric bootstrap on the point statistic from mean_relative_deviation.
    Returns (point_estimate, ci_low, ci_high) as percentages.
    """
    n = int(counts.sum())
    species = counts.index.to_numpy()
    p_hat = counts.to_numpy() / n

    resampled = rng.multinomial(n, p_hat, size=n_bootstrap)
    stats = np.empty(n_bootstrap)
    for i in range(n_bootstrap):
        stats[i] = mean_relative_deviation(pd.Series(resampled[i], index=species), expected)

    point = mean_relative_deviation(counts, expected)
    ci_low, ci_high = np.percentile(stats, [2.5, 97.5])
    return point, ci_low, ci_high


def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_CLASSIFIED_CSV
    if not os.path.exists(csv_path):
        raise SystemExit(f"No classification CSV at {csv_path}.")

    df = pd.read_csv(csv_path)
    classified = df[df["best_match"].notna()]
    rng = np.random.default_rng(SEED)

    # Winner-take-all: every classified read counts once, for whichever
    # species won its vote. Categories = all 10 species that can appear as
    # best_match, not just the 8 scored ones - see mean_relative_deviation's
    # docstring for why that distinction matters.
    wta_counts = classified["best_match"].value_counts()
    wta_point, wta_lo, wta_hi = bootstrap_ci(wta_counts, EXPECTED_ABUNDANCE, N_BOOTSTRAP, rng)
    print(f"Winner-take-all:   {wta_point:.1f}% [{wta_lo:.1f}%, {wta_hi:.1f}%]  (n={int(wta_counts.sum()):,} reads)")

    # Discard-ambiguous: same statistic, restricted to reads whose k-mers
    # hit exactly one species (n_species_hit == 1) - the subset Zymo's own
    # suggested remedy amounts to keeping.
    unambiguous = classified[classified["n_species_hit"] == 1]
    da_counts = unambiguous["best_match"].value_counts()
    da_point, da_lo, da_hi = bootstrap_ci(da_counts, EXPECTED_ABUNDANCE, N_BOOTSTRAP, rng)
    print(f"Discard-ambiguous: {da_point:.1f}% [{da_lo:.1f}%, {da_hi:.1f}%]  (n={int(da_counts.sum()):,} reads)")

    # Redistribution: point estimate only. A read's vote is split fractionally
    # across every species it hit, so outcomes are correlated within a read
    # rather than one independent categorical draw - the multinomial model
    # above doesn't hold, and it already loses by 12+ points anyway, so a CI
    # wouldn't change any conclusion.
    redistributed_path = Path(csv_path).with_name(Path(csv_path).stem + "_redistributed.csv")
    if redistributed_path.exists():
        redist_df = pd.read_csv(redistributed_path)
        redist_pct = dict(zip(redist_df["species"], redist_df["percentage"]))
        deviations = [abs(redist_pct.get(s, 0) - exp) / exp for s, exp in EXPECTED_ABUNDANCE.items()]
        redist_point = sum(deviations) / len(deviations) * 100
        print(f"Redistribution:    {redist_point:.1f}%  (point estimate only, see comment above)")
    else:
        print(f"Redistribution:    no redistributed CSV found at {redistributed_path}, skipping")


if __name__ == "__main__":
    main()