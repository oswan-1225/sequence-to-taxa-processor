
def parse_fasta(file_path: str) -> dict:
    """
    Parses a FASTA file and returns a dictionary of sequences
    
    Args:
        file_path (str): Path to the FASTA file

    Returns:
        dict: A dictionary of sequences
    """
    sequences = {}
    with open(file_path) as file:
        sequence_name = None
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence_name:
                    sequences[sequence_name] = ''.join(sequence)
                sequence_name = line[1:]  # Remove the '>' character
                sequence = []
            else:
                sequence.append(line)
        if sequence_name:
            sequences[sequence_name] = ''.join(sequence)
    return sequences

def extract_kmers(sequence: str, k: int) -> list:
    """
    Extracts k-mers from a given sequence.
    
    Args:
        sequence (str): The input sequence
        k (int): The length of the k-mers to extract

    Returns:
        list: A list of k-mers
    """
    if k <= 0:
        raise ValueError("k must be a positive integer")
    if k > len(sequence):
        return []
    
    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    return kmers

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
    for seq_name, seq in sequences.items():
        kmers = extract_kmers(seq, k)
        for kmer in kmers:
            kmer_index.setdefault(kmer, set()).add(seq_name)
    return kmer_index