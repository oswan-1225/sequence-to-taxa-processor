import argparse
from typing import Optional
from fasta_utils import parse_sequence_file, parse_fastq_qualities
from classifier_functions import classify_read_top_hit, load_index
from qc import filter_low_quality_reads
from tqdm import tqdm
import pandas as pd
from pathlib import Path
from database import insert_sample_results, create_database
from collections import defaultdict

def classify_file(index_path: str, reads_path: str, k: int, output_path: str, db_path: str = None,
                   source: str = None, min_quality: Optional[float] = None,
                   redistribute: bool = False) -> pd.DataFrame:
    """
    Classify every read in reads_path against a saved k-mer index, save the
    results to output_path (.csv), and optionally insert them into a
    classification database.

    Parameters:
        min_quality: if given, drop reads whose mean Phred quality score is
            below this threshold before classifying (see
            qc.filter_low_quality_reads). Requires a FASTQ input - raises
            ValueError if reads_path is FASTA (no quality scores present).
        redistribute: if True, also write a second CSV
            (`<output_path>_redistributed.csv`) with each read's k-mer
            votes split proportionally across every species it hit,
            instead of winner-take-all. Off by default: on real 16S
            amplicon data this measures worse than winner-take-all (see
            README "Validated accuracy"), so it isn't the recommended
            method - it's kept as an opt-in for data where it may still
            help (e.g. true shotgun sequencing, where ambiguous
            multi-species reads are rare rather than the norm).

    Returns:
        pd.DataFrame: one row per read (read_id, best_match, confidence) -
            the same data written to output_path.
    """
    kmer_index = load_index(index_path, k)
    print(f"Loaded index with {len(kmer_index)} k-mers")

    reads = parse_sequence_file(reads_path)
    print(f"Loaded {len(reads)} reads from {reads_path}")

    if min_quality is not None:
        with open(reads_path) as f:
            first_char = f.read(1)
        if first_char != '@':
            raise ValueError("--min-quality filtering requires a FASTQ input "
                              "(quality scores aren't present in FASTA files).")
        qualities = parse_fastq_qualities(reads_path)
        before = len(reads)
        reads = filter_low_quality_reads(reads, qualities, min_quality)
        print(f"Quality filter (min_quality={min_quality}): kept {len(reads)}/{before} reads")

    vote_accumulator = defaultdict(float) if redistribute else None
    results = []

    for read_name, read_seq in tqdm(reads.items(), desc="Classifying reads"):
        classification = classify_read_top_hit(read_seq, kmer_index, k)
        if redistribute:
            votes = classification['votes']
            totalvotes = sum(votes.values())
            if votes:
                for species, count in votes.items():
                    vote_accumulator[species] += count / totalvotes

        results.append({
            "read_id": read_name,
            "best_match": classification['best_match'],
            "confidence": classification['confidence']
        })
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_path, index=False)
    print(f"Classification results saved to {output_path}")

    if redistribute:
        total_classified = sum(vote_accumulator.values())
        if total_classified > 0:
            output_path_obj = Path(output_path)
            redistributed_path = output_path_obj.with_name(output_path_obj.stem + "_redistributed.csv")
            redistributed_df = pd.DataFrame([
                {
                    "species": species,
                    "estimated_reads": credit,
                    "percentage": credit / total_classified * 100
                }
                for species, credit in vote_accumulator.items()
            ]).sort_values("percentage", ascending=False)
            redistributed_df.to_csv(redistributed_path, index=False)
            print(f"Redistributed abundance estimates saved to {redistributed_path}")

    if db_path:
        create_database(db_path)
        sample_id = Path(reads_path).stem.removesuffix("_1").removesuffix("_2")
        insert_sample_results(db_path, sample_id, source or "unknown", results)
        print(f"Inserted {len(results)} classification rows into database at {db_path} (sample_id={sample_id})")

    return results_df


def main():
    parser = argparse.ArgumentParser(description="Classify sequencing reads against a pre-built k-mer reference index.")
    parser.add_argument("--index", required=True, help="Path to a k-mer index built by build_reference.py (.pkl)")
    parser.add_argument("--reads", required=True, help="Path to a FASTA or FASTQ file of reads to classify")
    parser.add_argument("--k", type=int, default=21, help="K-mer length - MUST match the value used to build the index")
    parser.add_argument("--output", required=True, help="Path to save classification results (.csv)")
    parser.add_argument("--db", required=False, help="Path to the SQLite database file")
    parser.add_argument("--source", required=False, help="Source of the sample (e.g., 'SRA', 'local') for database entry")
    parser.add_argument("--min-quality", type=float, required=False, default=None,
                         help="Minimum mean Phred quality score to keep a read (FASTQ input only). "
                              "Reads below this are dropped before classification.")
    parser.add_argument("--redistribute", action="store_true", default=False,
                         help="Also write a second CSV splitting each read's k-mer votes "
                              "proportionally across every species it hit, instead of "
                              "winner-take-all. Off by default - not recommended for 16S "
                              "amplicon data, see README 'Validated accuracy'.")

    args = parser.parse_args()

    classify_file(args.index, args.reads, args.k, args.output, args.db, args.source, args.min_quality,
                  args.redistribute)

if __name__ == "__main__":
    main()
