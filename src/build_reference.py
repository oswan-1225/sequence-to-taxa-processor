import argparse
import os
from fasta_utils import parse_sequence_file
from classifier_functions import build_kmer_index
import pickle

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GENOME_DIR = os.path.join(SCRIPT_DIR, "..", "data", "reference", "genomes") # Path to the directory containing the genome .fna files
K = 21 # Kmer length
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "..", "data", "reference", "kmer_index.pkl") # Outputs the k-mer index to a pickle file for later use

VALID_EXTENSIONS = {'.fna', '.fasta', '.fa', '.fastq', '.fq'}

def load_all_genomes(genome_dir: str) -> dict:
    """
    Loads every recognized sequence file (FASTA or FASTQ, by extension) in a
    given directory into a dictionary of {species name: sequence}, merging
    multi-sequence genomes (chromosomes + plasmid) into a one species-level key.
    """

    all_sequences = {}
    for filename in os.listdir(genome_dir):
        ext = os.path.splitext(filename)[1].lower()
        if ext in VALID_EXTENSIONS:
            species_name = os.path.splitext(filename)[0]
            file_path = os.path.join(genome_dir, filename)
            sequences = parse_sequence_file(file_path)
            all_sequences[species_name] = ''.join(sequences.values())
    return all_sequences

def main():
    parser = argparse.ArgumentParser(description="Build a k-mer reference index from a folder of reference genomes.")
    parser.add_argument("--genome_dir", required=True, help="Path to a folder of reference genome files (.fna/.fasta/.fa/.fastq/.fq)")
    parser.add_argument("--output", required=True, help="Where to save the resulting k-mer index (.pkl)")
    parser.add_argument("--k", type=int, default=21, help="K-mer length (default: 21)")

    args = parser.parse_args()

    genomes = load_all_genomes(args.genome_dir)
    print(f"Loaded {len(genomes)} genomes from {args.genome_dir}")
    for name, seq in genomes.items():
        print(f"  {name}: {len(seq)} bases")

    index = build_kmer_index(genomes, args.k)
    print(f"Built k-mer index with {len(index)} unique k-mers of length {args.k}")

    with open(args.output, 'wb') as f:
        pickle.dump(index, f)
    print(f"Saved index to {args.output}")

if __name__ == "__main__":
    main()