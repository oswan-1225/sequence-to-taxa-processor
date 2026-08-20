import argparse
import time
from typing import Optional
from classifier_functions import classify_read_top_hit
from fasta_utils import parse_sequence_file
import pickle


def run_full_benchmark(reads: dict, index_path: str, k: int, sample_size: Optional[int] = None) -> dict:
    """
    End-to-end benchmark: load index + classify reads, measuring total wall-clock time.

    Parameters:
        reads: dict of {read_id: sequence}
        index_path: path to the pickled k-mer index
        k: k-mer length (must match the index)
        sample_size: if given, only classify the first N reads

    Returns:
        dict summarizing load time, classification time, and total time
    """
    overall_start = time.time()

    load_start = time.time()
    with open(index_path, "rb") as f:
        kmer_index = pickle.load(f)
    load_elapsed = time.time() - load_start

    read_items = list(reads.items())
    if sample_size is not None:
        read_items = read_items[:sample_size]

    classify_start = time.time()
    results = []
    for read_id, sequence in read_items:
        classification = classify_read_top_hit(sequence, kmer_index, k)
        results.append({"read_id": read_id, **classification})
    classify_elapsed = time.time() - classify_start

    total_elapsed = time.time() - overall_start
    n_reads = len(read_items)

    print(f"Index load time:        {load_elapsed:.2f}s ({len(kmer_index):,} k-mers)")
    print(f"Classification time:    {classify_elapsed:.2f}s ({n_reads:,} reads, "
          f"{n_reads/classify_elapsed:.1f} reads/sec)")
    print(f"Total end-to-end time:  {total_elapsed:.2f}s")

    return {
        "n_reads": n_reads,
        "load_seconds": load_elapsed,
        "classify_seconds": classify_elapsed,
        "total_seconds": total_elapsed,
        "reads_per_second": n_reads / classify_elapsed,
        "results": results
    }


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark index load + classification throughput end to end."
    )
    parser.add_argument("--index", required=True,
                        help="Path to a k-mer index built by build_reference.py (.pkl)")
    parser.add_argument("--reads", required=True,
                        help="Path to a FASTA or FASTQ file of reads to classify")
    parser.add_argument("--k", type=int, default=21,
                        help="K-mer length - MUST match the value used to build the index")
    parser.add_argument("--sample-size", type=int, default=None,
                        help="Only classify the first N reads (default: all of them)")

    args = parser.parse_args()

    reads = parse_sequence_file(args.reads)
    run_full_benchmark(reads, args.index, args.k, args.sample_size)


if __name__ == "__main__":
    main()
