
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

sequences = parse_fasta("data/reference/genomes/Escherichia_coli_complete_genome.fna")
print(len(sequences))
for name, seq in sequences.items():
    print(name, len(seq))