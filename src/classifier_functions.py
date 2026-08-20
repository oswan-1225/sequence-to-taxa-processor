from fasta_utils import parse_fasta
from fasta_utils import extract_kmers
from tqdm import tqdm

def build_kmer_index(sequences: dict, k: int) -> dict:
    """
    Builds an index of k-mers from a dict of sequences
    
    Args:
        sequences (dict): A dictionary of sequences
        k (int): The length of the k-mers to extract
        
    Returns:
        dict: A dictionary keys: k-mers, values: list of sequence names containing the k-mer
    """
    kmer_index = {}
    for seq_name, seq in tqdm(sequences.items(), desc="Indexing genomes"):
        kmers = extract_kmers(seq, k)
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
    read_kmers = extract_kmers(read_sequence, k)
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
    total_kmers = len(extract_kmers(read_sequence, k))

    if not classif or total_kmers == 0:
        return {'best_match': None, 'confidence': 0.0, 'votes': classif}

    top_count = max(classif.values())
    winners = [species for species, count in classif.items() if count == top_count]

    if len(winners) > 1:
        return {'best_match': None, 'confidence': 0.0, 'votes': classif}

    return {'best_match': winners[0], 'confidence': top_count / total_kmers, 'votes': classif}