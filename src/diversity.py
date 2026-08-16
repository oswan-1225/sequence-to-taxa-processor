import numpy as np
import pandas as pd


def species_richness(sample_counts: list | pd.Series, threshold: float = 0.0) -> int:
    """
    Number of species present (count > threshold) in a single sample.

    Parameters:
        sample_counts: array/Series of per-species counts for one sample
        threshold: minimum count for a species to be considered "present"

    Returns:
        int: number of species with count > threshold
    """
    return int(np.sum(np.array(sample_counts) > threshold))


def shannon_diversity(sample_counts: list | pd.Series) -> float:
    """
    Shannon diversity index (H = -sum(p_i * ln(p_i))) for a single sample.

    Parameters:
        sample_counts: array/Series of per-species counts for one sample

    Returns:
        float: Shannon diversity index. 0.0 for an empty/all-zero sample.
    """
    counts = np.array(sample_counts)
    if np.all(counts == 0):
        return 0.0

    total = np.sum(counts)
    if total == 0:
        return 0.0
  
    proportions = counts / total
    nonzero = proportions > 0
    proportions = proportions[nonzero]
    
    log_proportions = np.log(proportions)
    return -np.sum(proportions * log_proportions)


def diversity_by_sample(abundance_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-sample diversity metrics from get_abundance()'s long/tidy output.

    Parameters:
        abundance_df: DataFrame from database.get_abundance(), with columns
            ['sample_id', 'best_match', 'count', 'total', 'percent'] - one
            row per (sample, species).

    Returns:
        pd.DataFrame: one row per sample_id with its diversity metric(s).
    """

    diversity_metrics = []
    for sample_id, group in abundance_df.groupby("sample_id"):
        counts = group['count']
        richness = species_richness(counts)
        shannon = shannon_diversity(counts)
        diversity_metrics.append({
            'sample_id': sample_id,
            'species_richness': richness,
            'shannon_diversity': shannon
        })

    return pd.DataFrame(diversity_metrics)
