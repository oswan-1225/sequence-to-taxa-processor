"""
Downloads the real Illumina reads used to validate this tool against the
ZymoBIOMICS D6300 mock community (SRA accession SRR10391187). Lives in
docs/ alongside plot_zymobiomics_validation.py - dataset-specific, not
part of the general-purpose pipeline.

Requires the SRA Toolkit (prefetch + fasterq-dump) on PATH.
"""
import os
import shutil
import subprocess
import sys

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ACCESSION = "SRR10391187"
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "sra_reads")
SRA_CACHE_DIR = os.path.join(OUTPUT_DIR, "_sra_cache")


def download_reads(accession: str = ACCESSION, output_dir: str = OUTPUT_DIR) -> None:
    for tool in ("prefetch", "fasterq-dump"):
        if shutil.which(tool) is None:
            raise SystemExit(
                f"'{tool}' not found on PATH. Install the SRA Toolkit "
                "(https://github.com/ncbi/sra-tools/wiki/HowTo:-Binary-Installation) "
                "and re-run this script."
            )

    fastq_1 = os.path.join(output_dir, f"{accession}_1.fastq")
    fastq_2 = os.path.join(output_dir, f"{accession}_2.fastq")
    if os.path.exists(fastq_1) and os.path.exists(fastq_2):
        print(f"{fastq_1} and {fastq_2} already exist, skipping download.")
        return

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(SRA_CACHE_DIR, exist_ok=True)

    print(f"Prefetching {accession}...")
    subprocess.run(["prefetch", accession, "-O", SRA_CACHE_DIR], check=True)

    sra_file = os.path.join(SRA_CACHE_DIR, accession, f"{accession}.sra")
    print(f"Extracting FASTQ from {sra_file}...")
    subprocess.run(["fasterq-dump", sra_file, "-O", output_dir, "--split-files"], check=True)

    print(f"Saved {fastq_1} and {fastq_2}")


def main():
    try:
        download_reads()
    except subprocess.CalledProcessError as e:
        sys.exit(f"'{e.cmd[0]}' failed with exit code {e.returncode}")


if __name__ == "__main__":
    main()
