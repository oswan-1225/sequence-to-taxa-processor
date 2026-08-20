
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
            f"({gc:.1f}%, z={z:+.2f} vs. set mean, threshold={threshold}). "
            f"This describes the reference genome, not a problem detected in "
            f"your data. It affects abundance only if the library prep "
            f"introduced GC-dependent coverage bias, which applies to "
            f"whole-genome shotgun sequencing but not to amplicon data, where "
            f"every read comes from one locus. Check your library type before "
            f"treating this species' estimate as less reliable."
        )
    return warnings

def mean_phred_quality(quality_string: str) -> float:
    """
    Calculate the mean Phred+33 quality score of a FASTQ quality string.
    Phred+33 (Illumina 1.8+, the standard encoding for all current SRA data):
    Q = ord(char) - 33, where Q estimates -10*log10(base call error
    probability) - e.g. Q20 = 1% error rate, Q30 = 0.1% error rate.

    Parameters:
    quality_string (str): raw ASCII quality line from a FASTQ record.

    Returns:
    float: mean Phred quality score across all bases. Empty string returns 0.0.
    """
    if not quality_string:
        return 0.0
    scores = [ord(char) - 33 for char in quality_string]
    return sum(scores) / len(scores)

def per_read_mean_quality(qualities: dict) -> dict:
    """
    Calls mean_phred_quality() per read and returns {read_id: mean quality}.

    Parameters:
    qualities (dict): {read_id: quality_string}, e.g. from
        fasta_utils.parse_fastq_qualities().

    Returns:
    dict: {read_id: mean Phred quality score}.
    """
    return {read_id: mean_phred_quality(q) for read_id, q in qualities.items()}

def filter_low_quality_reads(reads: dict, qualities: dict, min_quality: float = 20.0) -> dict:
    """
    Drop reads whose mean Phred quality score is below min_quality.
    Default of 20.0 (Q20, ~99% per-base call accuracy) is the standard
    baseline threshold used across NGS QC tools (e.g. FastQC/Trimmomatic).

    Parameters:
    reads (dict): {read_id: sequence}.
    qualities (dict): {read_id: quality_string}, same read_ids as reads.
    min_quality (float): minimum mean Phred quality score required to keep a read.

    Returns:
    dict: {read_id: sequence}, filtered to only reads meeting min_quality.
    """
    mean_qualities = per_read_mean_quality(qualities)
    return {
        read_id: seq for read_id, seq in reads.items()
        if mean_qualities.get(read_id, 0.0) >= min_quality
    }