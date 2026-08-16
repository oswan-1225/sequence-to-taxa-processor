
def gc_content(sequence: str) -> float:
    """
    Calculate the GC content of a DNA sequence.

    Parameters:
    sequence (str): The DNA sequence.

    Returns:
    float: The GC content as a percentage.
    """
    if not sequence:
        return 0.0
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100

def gc_zscores(gc_by_species: dict) -> dict:
    """
    Calculate how many standard deviations each species' GC content is from
    the mean GC content across all species in the given reference set.

    Parameters:
    gc_by_species (dict): {species name: GC content percentage}, one value
        per species - e.g. gc_content() called once per genome.

    Returns:
    dict: {species name: z-score}. 0.0 means "average for this set"; larger
        magnitude (either sign) means further from the mean.
    """
    values = list(gc_by_species.values())
    if len(values) < 2:
        return {species: 0.0 for species in gc_by_species}

    mean = sum(values) / len(values)
    variance = sum((x - mean) ** 2 for x in values) / (len(values) - 1)  # sample std (ddof=1)
    std_dev = variance ** 0.5

    if std_dev == 0:
        return {species: 0.0 for species in gc_by_species}

    return {species: (gc - mean) / std_dev for species, gc in gc_by_species.items()}

def flag_gc_outliers(gc_by_species: dict, threshold: float = 1.5) -> dict:
    """
    Flag species whose GC content is more than 'threshold' std away from mean GC content across all species in the given reference set.

    Parameters:
    gc_by_species (dict): {species name: GC content percentage}, one value
        per species - e.g. gc_content() called once per genome.
    threshold (float): Number of standard deviations away from the mean to consider as an outlier.

    Returns:
    dict: {species name: True/False}, where True indicates the species is an outlier.
    """
    zscores = gc_zscores(gc_by_species)
    return {species: abs(z) > threshold for species, z in zscores.items()}

def gc_outlier_warnings(gc_by_species: dict, threshold: float = 1.5) -> dict:
    """
    Generate a warning message for each species whose GC content is a
    statistical outlier within the given reference set.

    Parameters:
    gc_by_species (dict): {species name: GC content percentage}, one value
        per species - e.g. gc_content() called once per genome.
    threshold (float): Number of standard deviations from the mean to
        consider as an outlier.

    Returns:
    dict: {species name: warning message}, containing only species flagged
        as outliers - non-outlier species are omitted entirely, so a caller
        can do e.g. df['col'].map(warnings) and get NaN for the rest.
    """
    zscores = gc_zscores(gc_by_species)
    flags = flag_gc_outliers(gc_by_species, threshold)

    warnings = {}
    for species, is_outlier in flags.items():
        if not is_outlier:
            continue
        z = zscores[species]
        gc = gc_by_species[species]
        direction = "high" if z > 0 else "low"
        warnings[species] = (
            f"GC content is unusually {direction} for this reference set "
            f"({gc:.1f}%, z={z:+.2f} vs. set mean, threshold={threshold}) - "
            f"sequencing/library-prep GC bias is known to under- or "
            f"over-represent extreme-GC genomes, so this species' abundance "
            f"estimate may be less reliable than the others."
        )
    return warnings

def per_read_gc_content(reads: dict) -> dict:
    """
    Calls gc_content() per read and returns a dict of {read_id: gc_content}. for histogramming.
    
    Parameters:
    reads (dict): {read_id: sequence} for all reads in a FASTA/FASTQ file.

    Returns:
    dict: {read_id: GC content percentage}.
    """
    return {read_id: gc_content(seq) for read_id, seq in reads.items()}

def reads_in_gc_range(reads: dict, target_gc: float, n_sds: float) -> int:
    """
    Count how many reads have GC content consistent with having been sampled
    from a genome with GC content target_gc, allowing for binomial sampling
    noise: a read of length n is a small random sample of its source genome,
    so its own GC% naturally scatters around the true genome GC% with
    SD = sqrt(target_gc * (100 - target_gc) / n). Each read uses its own
    length, so the tolerance window is wider for shorter reads and narrower
    for longer ones - not one fixed window shared by every read.

    Parameters:
    reads (dict): {read_id: sequence} for all reads in a FASTA/FASTQ file.
    target_gc (float): Target genome GC content percentage to compare against.
    n_sds (float): Number of standard deviations from target_gc to allow.

    Returns:
    int: Count of reads within the specified GC content range.
    """
    count = 0
    for sequence in reads.values():
        n = len(sequence)
        if n == 0:
            continue

        sd = (target_gc * (100 - target_gc) / n) ** 0.5
        lower_bound = target_gc - n_sds * sd
        upper_bound = target_gc + n_sds * sd

        if lower_bound <= gc_content(sequence) <= upper_bound:
            count += 1

    return count