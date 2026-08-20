import pickle

from fasta_utils import parse_fasta
from fasta_utils import extract_kmers, extract_canonical_kmers
from tqdm import tqdm

# Bumped whenever a saved index's contents stop being interchangeable with
# older ones. Version 1 was a bare {kmer: set(species)} dict of FORWARD-strand
# k-mers only, with no record of k or strandedness. Version 2 wraps that dict
# in a header and canonicalizes every k-mer.
#
# The version exists because reusing a v1 index against canonical queries does
# not fail, nearly every lookup simply misses, and you get a complete,
# plausible-looking CSV of nulls. An index is an opaque 2.2GB binary that takes
# minutes to build, so "just rebuild it to be safe" is not a cheap habit, which
# means the file has to be able to say what it is.
INDEX_FORMAT_VERSION = 2


def save_index(path: str, kmer_index: dict, k: int) -> None:
    """
    Write a k-mer index to disk with the metadata needed to validate it later.

    Args:
        path (str): Where to write the .pkl
        kmer_index (dict): The index from build_kmer_index
        k (int): The k-mer length it was built with
    """
    payload = {
        "format_version": INDEX_FORMAT_VERSION,
        "k": k,
        "canonical": True,
        "index": kmer_index,
    }
    with open(path, 'wb') as f:
        pickle.dump(payload, f)


def load_index(path: str, k: int) -> dict:
    """
    Read a k-mer index from disk, refusing anything that would silently
    produce wrong results.

    k is required rather than optional so the check cannot be skipped by
    accident. A k mismatch has always been a silent failure in this tool: the
    CLIs take --k independently of the index, and nothing verified the two
    agreed. A read tokenized at k=25 against an index built at k=21 matches
    almost nothing, and the run completes normally.

    Args:
        path (str): Path to a .pkl written by save_index
        k (int): The k-mer length the caller intends to use

    Returns:
        dict: The {kmer: set(species)} index itself, unwrapped

    Raises:
        ValueError: if the file predates the format, was written by a newer
            version, or was built with a different k
    """
    with open(path, 'rb') as f:
        payload = pickle.load(f)

    # A v1 index is a bare dict keyed on DNA strings, so "format_version"
    # cannot collide with a real key - only the letters ACGT appear in one.
    if not isinstance(payload, dict) or "format_version" not in payload:
        raise ValueError(
            f"{path} is an old-format k-mer index (forward-strand only, "
            f"pre-canonical). Its results would be wrong rather than merely "
            f"stale, so it cannot be reused. Rebuild it with:\n"
            f"  python src/build_reference.py --genome_dir <genomes> "
            f"--output {path} --k {k}"
        )

    if payload["format_version"] != INDEX_FORMAT_VERSION:
        raise ValueError(
            f"{path} was written in index format v{payload['format_version']}, "
            f"but this code reads v{INDEX_FORMAT_VERSION}. Rebuild the index "
            f"with src/build_reference.py."
        )

    if payload["k"] != k:
        raise ValueError(
            f"{path} was built with k={payload['k']}, but k={k} was requested. "
            f"Reads tokenized at a different k than the index match almost "
            f"nothing. Pass --k {payload['k']}, or rebuild the index at k={k}."
        )

    return payload["index"]


def build_kmer_index(sequences: dict, k: int) -> dict:
    """
    Builds an index of k-mers from a dict of sequences
    
    Args:
        sequences (dict): A dictionary of sequences
        k (int): The length of the k-mers to extract
        
    Returns:
        dict: A dictionary keys: k-mers, values: set of sequence names containing the k-mer
    """
    kmer_index = {}
    for seq_name, seq in tqdm(sequences.items(), desc="Indexing genomes"):
        kmers = extract_canonical_kmers(seq, k)
        for kmer in kmers:
            kmer_index.setdefault(kmer, set()).add(seq_name)
    return kmer_index

def classify_read(read_sequence: str, kmer_index: dict, k: int) -> dict:
    """
    Classifies a read based on the k-mer index.
    
    Args:
        read_sequence (str): The sequence of the read
        kmer_index (dict): The k-mer index
        k (int): The length of the k-mers to extract
        
    Returns:
        dict: A dictionary of sequences names and their k-mer counts in the read
    """
    read_kmers = extract_canonical_kmers(read_sequence, k)
    classification = {}
    for kmer in read_kmers:
        if kmer in kmer_index:
            for seq_name in kmer_index[kmer]:
                classification[seq_name] = classification.get(seq_name, 0) + 1
    return classification

def classify_read_top_hit(read_sequence: str, kmer_index: dict, k: int) -> dict:
    """
    Classifies a read and returns its single best matching sequence, if it has one.

    A read is assigned to a species only when that species holds the top k-mer
    vote count outright. When two or more species tie for the lead the read is
    left unassigned (best_match None) rather than awarded to one of them: the
    classifier genuinely cannot tell them apart, and any tie-break would either
    bias a species systematically or depend on set/dict iteration order, which
    varies between processes and makes runs non-reproducible. Measured on
    SRR10391187 (k=21, 10-species index), ties affect ~0.15% of reads.

    Args:
        read_sequence (str): The sequence of the read
        kmer_index (dict): The k-mer index
        k (int): The length of the k-mers to extract

    Returns:
        dict: 'best_match' (species name, or None if the read matched nothing
            or tied for the lead), 'confidence' (best_match's share of the
            read's k-mers, 0.0 whenever best_match is None) and 'votes' (the
            full {species: kmer_hit_count} tally, populated even for a tie so
            callers can see what the read actually hit)
    """
    classif = classify_read(read_sequence, kmer_index, k)
    # Deliberately extract_kmers, not extract_canonical_kmers. This is only a
    # denominator - the number of k-mer POSITIONS in the read - and
    # canonicalizing rewrites k-mers without adding or dropping positions, so
    # the count is identical either way. Canonicalizing here would double this
    # function's reverse-complement work for no change in result.
    total_kmers = len(extract_kmers(read_sequence, k))

    if not classif or total_kmers == 0:
        return {'best_match': None, 'confidence': 0.0, 'votes': classif}

    top_count = max(classif.values())
    winners = [species for species, count in classif.items() if count == top_count]

    if len(winners) > 1:
        return {'best_match': None, 'confidence': 0.0, 'votes': classif}

    return {'best_match': winners[0], 'confidence': top_count / total_kmers, 'votes': classif}