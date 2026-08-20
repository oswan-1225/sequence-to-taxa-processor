import time

# Table for making complementary k-mers
_COMPLEMENT = str.maketrans({
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
    'a': 't',
    't': 'a',
    'c': 'g',
    'g': 'c',
})

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

def reverse_complement(sequence: str) -> str:
    """ Returns the reverse complement of a DNA sequence """
    complemented = sequence.translate(_COMPLEMENT)

    reversed_sequence = complemented[::-1]

    return reversed_sequence

def canonical_kmer(kmer: str) -> str:
    """ Returns whichever of kmer and its reverse complement sorts first """
    rc = reverse_complement(kmer)
    return kmer if kmer <= rc else rc 

def extract_canonical_kmers(sequence: str, k: int) -> list:
    """
    Extracts every k-mer in a sequence, in canonical form.

    A sequencer reports whichever strand of a double-stranded fragment it
    happened to read, and does not say which. The same physical DNA can
    therefore arrive written either as a k-mer or as that k-mer's reverse
    complement. Both forms have to collapse onto one dictionary key or the
    index cannot match them, which is what canonical_kmer does.

    This must be applied to BOTH sides - the reference genomes in
    build_kmer_index and the reads in classify_read. Canonicalizing only one
    side would make every lookup miss.

    Delegates to extract_kmers rather than reimplementing the slicing, so
    there is exactly one definition of "the k-mers of a sequence" in the
    codebase, and its k <= 0 and k > len(sequence) handling is inherited
    rather than duplicated.

    Args:
        sequence (str): The input sequence
        k (int): The length of the k-mers to extract

    Returns:
        list: The canonical k-mers, in sequence order. Same length as
            extract_kmers(sequence, k) - canonicalizing rewrites k-mers, it
            never adds or drops positions.
    """
    return [canonical_kmer(kmer) for kmer in extract_kmers(sequence, k)]



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


def parse_fastq_qualities(filepath: str) -> dict:
    """
    Parses a FASTQ file's quality lines into a dict of {read_id: quality_string}.
    Companion to parse_fastq(), which parses the same file's sequences but
    discards quality - kept as a separate function (rather than changing
    parse_fastq's return shape) so every existing caller of parse_fastq /
    parse_sequence_file is unaffected.

    Parameters:
        filepath: path to a .fastq file

    Returns:
        dict: {read_id: quality_string} - the raw ASCII quality line,
            Phred+33 encoded (see qc.mean_phred_quality for decoding).
    """
    qualities = {}
    with open(filepath) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break # End of file
            f.readline() # Skip the sequence line
            f.readline() # Skip the '+' line
            quality = f.readline().strip()

            read_id = header[1:] # remove '@' from header
            qualities[read_id] = quality
    return qualities


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

