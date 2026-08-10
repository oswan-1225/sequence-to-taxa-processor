import argparse
import pickle
from fasta_utils import parse_fasta
from classifier_functions import classify_read_top_hit

def main():
    parser = argparse.ArgumentParser(description="Classify sequencing reads against a pre-built k-mer reference index.")
    parser.add_argument("--index", required=True, help="Path to a k-mer index built by build_reference.py (.pkl)")
    parser.add_argument("--reads", required=True, help="Path to a FASTA file of reads to classify")
    parser.add_argument("--k", type=int, default=21, help="K-mer length - MUST match the value used to build the index")
    parser.add_argument("--output", required=True, help="Path to save classification results (.csv)")

    args = parser.parse_args()

    with open(args.index, 'rb') as f:
        kmer_index = pickle.load(f)
    print(f"Loaded index with {len(kmer_index)} k-mers")

    reads = parse_fasta(args.reads)
    print(f"Loaded {len(reads)} reads from {args.reads}")


    results = []
    for read_name, read_seq in reads.items():
        classification = classify_read_top_hit(read_seq, kmer_index, args.k)
        results.append((read_name, classification['best_match'], classification['confidence']))   
if __name__ == "__main__":
    main()