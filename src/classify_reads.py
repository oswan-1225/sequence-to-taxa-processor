import argparse
import pickle
from fasta_utils import parse_sequence_file
from classifier_functions import classify_read_top_hit
from tqdm import tqdm
import pandas as pd
from pathlib import Path
from database import insert_sample_results, create_database

def classify_file(index_path: str, reads_path: str, k: int, output_path: str, db_path: str = None, source: str = None) -> pd.DataFrame:
    """
    Classify every read in reads_path against a saved k-mer index, save the
    results to output_path (.csv), and optionally insert them into a
    classification database.

    Returns:
        pd.DataFrame: one row per read (read_id, best_match, confidence) -
            the same data written to output_path.
    """
    with open(index_path, 'rb') as f:
        kmer_index = pickle.load(f)
    print(f"Loaded index with {len(kmer_index)} k-mers")

    reads = parse_sequence_file(reads_path)
    print(f"Loaded {len(reads)} reads from {reads_path}")

    results = []
    for read_name, read_seq in tqdm(reads.items(), desc="Classifying reads"):
        classification = classify_read_top_hit(read_seq, kmer_index, k)
        results.append({
            "read_id": read_name,
            "best_match": classification['best_match'],
            "confidence": classification['confidence']
        })
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_path, index=False)
    print(f"Classification results saved to {output_path}")

    if db_path:
        create_database(db_path)
        sample_id = Path(reads_path).stem.removesuffix("_1").removesuffix("_2")
        insert_sample_results(db_path, sample_id, source or "unknown", results)
        print(f"Inserted {len(results)} classification rows into database at {db_path} (sample_id={sample_id})")

    return results_df


def main():
    parser = argparse.ArgumentParser(description="Classify sequencing reads against a pre-built k-mer reference index.")
    parser.add_argument("--index", required=True, help="Path to a k-mer index built by build_reference.py (.pkl)")
    parser.add_argument("--reads", required=True, help="Path to a FASTA file of reads to classify")
    parser.add_argument("--k", type=int, default=21, help="K-mer length - MUST match the value used to build the index")
    parser.add_argument("--output", required=True, help="Path to save classification results (.csv)")
    parser.add_argument("--db", required=False, help="Path to the SQLite database file")
    parser.add_argument("--source", required=False, help="Source of the sample (e.g., 'SRA', 'local') for database entry")

    args = parser.parse_args()

    classify_file(args.index, args.reads, args.k, args.output, args.db, args.source)

if __name__ == "__main__":
    main()
