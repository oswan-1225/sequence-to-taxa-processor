"""
One-off investigation (commit 3): is Pseudomonas aeruginosa's low classified
abundance (4.9% vs. ~12% expected) explained by raw-read undersequencing
consistent with its GC-atypical genome (66%+ GC, a real outlier among the
8 ZymoBIOMICS reference genomes - see src/qc.py's gc_zscores), or does it
point at something else (e.g. the classifier's own vote-sharing logic)?

Dataset-specific - the species name, file paths, and interpretation below
are all specific to this project's ZymoBIOMICS validation run, which is why
this lives here and not in src/qc.py (general-purpose library code takes
target_gc/n_sds as arguments; it doesn't know about "Pseudomonas" or
"SRR10391187").
"""
import os
import sys

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "src"))

from build_reference import load_all_genomes
from fasta_utils import parse_sequence_file
from qc import gc_content, reads_in_gc_range

GENOME_DIR = os.path.join(PROJECT_ROOT, "data", "reference", "Genomes")
READS_PATH = os.path.join(PROJECT_ROOT, "data", "raw", "sra_reads", "SRR10391187_1.fastq")
TARGET_SPECIES = "Pseudomonas_aeruginosa_complete_genome"

# ~95% tolerance under a normal approximation to the binomial sampling
# distribution - a read is counted as "consistent with target_gc" if its
# own GC% falls within +/-2 of its own expected sampling SD.
N_SDS = 2.0

# 8 species in the ZymoBIOMICS community, assumed evenly represented in the
# raw reads themselves (not the same claim as "evenly classified" - that's
# exactly what this script is checking).
EXPECTED_FRACTION = 1 / 8


def main():
    print("Loading reference genomes...")
    genomes = load_all_genomes(GENOME_DIR)
    gc_by_species = {species: gc_content(seq) for species, seq in genomes.items()}

    if TARGET_SPECIES not in gc_by_species:
        raise KeyError(f"{TARGET_SPECIES!r} not found in {GENOME_DIR}. "
                        f"Available species: {sorted(gc_by_species)}")
    target_gc = gc_by_species[TARGET_SPECIES]

    print(f"\nReference genome GC%, sorted low to high "
          f"(checking whether any other species' GC is close enough to "
          f"{TARGET_SPECIES}'s to land in the same window and confound the count):")
    for species, gc in sorted(gc_by_species.items(), key=lambda kv: kv[1]):
        marker = "  <-- target" if species == TARGET_SPECIES else ""
        print(f"  {gc:6.2f}%  {species}{marker}")

    print(f"\nLoading reads from {READS_PATH} ...")
    reads = parse_sequence_file(READS_PATH)
    total_reads = len(reads)
    print(f"Loaded {total_reads} reads.")

    print(f"\nCounting reads within {N_SDS} SDs of target_gc={target_gc:.2f}% "
          f"(each read's SD computed from its own length, per qc.reads_in_gc_range)...")
    observed_count = reads_in_gc_range(reads, target_gc, N_SDS)
    observed_fraction = observed_count / total_reads

    print(f"\nObserved: {observed_count}/{total_reads} reads "
          f"({observed_fraction:.1%}) have GC% consistent with a "
          f"{target_gc:.2f}%-GC source, within {N_SDS} SDs of sampling noise.")
    print(f"Expected if all 8 species were evenly represented in raw reads: "
          f"{EXPECTED_FRACTION:.1%}")
    print(f"Ratio (observed / expected): {observed_fraction / EXPECTED_FRACTION:.2f}x")


if __name__ == "__main__":
    main()
