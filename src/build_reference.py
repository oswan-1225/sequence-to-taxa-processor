import os
from fasta_utils import parse_fasta
from classifier import build_kmer_index

GENOME_DIR = 'data/reference/genomes'
K = 21 # Kmer length
OUTPUT_FILE = 'data/reference/kmer_index.pkl'

def load_all_genomes(genome_dir: str) -> dict:
    """
    Loads every .fna file in a given directory into a dictiojnary of 
    {species name: sequence}, merging multi-sequence genomes
    (chromosomes + plasmid} into a one species level key.)
    """

    all_sequences = {}
    for filename in os.listdir(genome_dir):
        if filename.endswith('fna'):
            species_name = os.path.splitext(filename)[0]
            file_path = os.path.join(genome_dir, filename)
            sequences = parse_fasta(file_path)
            all_sequences[species_name] = ''.join(sequences.values())
    return all_sequences

if __name__ == "__main__":
    genomes = load_all_genomes(GENOME_DIR)
    print(f"Loaded {len(genomes)} genomes from {GENOME_DIR}")
    for name, seq in genomes.items():
        print(f"{name}: {len(seq)} bases")

    index = build_kmer_index(genomes, K)
    print(f"Built k-mer index with {len(index)} unique k-mers of length {K}")