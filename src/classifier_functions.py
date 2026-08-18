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
    Classifies a read and returns the single best matching sequence based on the k-mer index with a confidence score.
    
    Args:
        read_sequence (str): The sequence of the read
        kmer_index (dict): The k-mer index
        k (int): The length of the k-mers to extract
    
    Returns:
        dict: A dictionary with the best matching sequence name and its confidence score
    """
    classif = classify_read(read_sequence, kmer_index, k)
    total_kmers = len(extract_kmers(read_sequence, k))

    if not classif or total_kmers == 0:
        return {'best_match': None, 'confidence': 0.0, 'votes': classif}
    best_match = max(classif.items(), key=lambda item: item[1])[0]
    return {'best_match': best_match, 'confidence': classif[best_match] / total_kmers, 'votes': classif}