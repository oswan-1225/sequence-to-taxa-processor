
import time


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

def parse_fastq(filepath: str) -> dict:
    """
    Parses a FASTQ file into a dict of {header: sequence}.
    Quality scores are read but discarded (not necessary for classification).

    Parameteres:
        filepath: path to a .fastq file

    Returns:
        dict: {read_header: read_sequence}
    """
    sequences = {}
    with open(filepath) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break # End of file
            sequence = f.readline().strip()
            f.readline() # Skip the '+' line
            f.readline() # Skip the quality score line

            read_id = header[1:] # rmove '@' from header
            sequences[read_id] = sequence
    return sequences

    
def parse_sequence_file(filepath: str) -> dict:
    """
    Parse a sequence file (FASTA or FASTQ)
    Prints timing and throughput info (reads/second) for benchmarking
    
    Returns:
        dict: {read_id: sequence}, regardless of format
    """
    start_time = time.time()

    with open(filepath) as f:
        first_char = f.read(1)

    if first_char == '>':
        sequences = parse_fasta(filepath)
    elif first_char == '@':
        sequences = parse_fastq(filepath)
    else:
        raise ValueError(f"File format not recognized. (starts with {first_char})")

    elapsed_time = time.time() - start_time
    rps = len(sequences) / elapsed_time if elapsed_time > 0 else float('inf')
    print(f"Parsed {len(sequences)} sequences in {elapsed_time:.2f} seconds, Throughput: {rps:.2f} reads/second")

    return sequences