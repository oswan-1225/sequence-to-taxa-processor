import argparse
import pickle
from fasta_utils import parse_fasta
from classifier_functions import classify_read_top_hit
from tqdm import tqdm
import pandas as pd
from pathlib import Path
from database import insert_sample_results, create_database

def main():
    parser = argparse.ArgumentParser(description="Classify sequencing reads against a pre-built k-mer reference index.")
    parser.add_argument("--index", required=True, help="Path to a k-mer index built by build_reference.py (.pkl)")
    parser.add_argument("--reads", required=True, help="Path to a FASTA file of reads to classify")
    parser.add_argument("--k", type=int, default=21, help="K-mer length - MUST match the value used to build the index")
    parser.add_argument("--output", required=True, help="Path to save classification results (.csv)")
    parser.add_argument("--db", required=False, help="Path to the SQLite database file")
    parser.add_argument("--source", required=False, help="Source of the sample (e.g., 'SRA', 'local') for database entry")

    args = parser.parse_args()

    with open(args.index, 'rb') as f:
        kmer_index = pickle.load(f)
    print(f"Loaded index with {len(kmer_index)} k-mers")

    reads = parse_fasta(args.reads)
    print(f"Loaded {len(reads)} reads from {args.reads}")

    results = []
    for read_name, read_seq in tqdm(reads.items(), desc="Classifying reads"):
        classification = classify_read_top_hit(read_seq, kmer_index, args.k)
        results.append({
            "read_id": read_name,
            "best_match": classification['best_match'],
            "confidence": classification['confidence']
        })
    results_df = pd.DataFrame(results)
    results_df.to_csv(args.output, index=False)
    print(f"Classification results saved to {args.output}")

    if args.db:
        create_database(args.db)
        sample_id = Path(args.reads).stem.removesuffix("_1").removesuffix("_2")
        insert_sample_results(args.db, sample_id, args.source or "unknown", results)
        print(f"Inserted {len(results)} classification rows into database at {args.db} (sample_id={sample_id})")

if __name__ == "__main__":
    main()
