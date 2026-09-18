import os
from typing import Optional

import click

from build_reference import load_all_genomes
from classifier_functions import build_kmer_index, save_index
from classify_reads import classify_file, derive_sample_id
from database import get_abundance, get_classification_totals
from diversity import diversity_by_sample
from qc import gc_content, gc_outlier_warnings
from visualization import plot_sample_summary


def run_pipeline(genome_dir: str, reads: str, output_dir: str, k: int = 21,
                  index: Optional[str] = None, source: Optional[str] = None, top_n: int = 10,
                  skip_plot: bool = False, min_quality: Optional[float] = None,
                  redistribute: bool = False, db: Optional[str] = None) -> dict:
    """
    Run the full build -> classify -> diversity -> plot pipeline for a single
    read file, writing every output into output_dir.

    Fails fast: each stage runs only if the previous one succeeded, and a
    failure is re-raised as a RuntimeError naming which stage failed.

    Parameters:
        genome_dir: folder of reference genomes. Required unless index is
            given. If provided (with or without index), it's also used to
            check reference genomes for GC-content outliers, regardless of
            whether it's needed to build the index - loading genome
            sequences is cheap relative to building the k-mer index, so this
            check runs even when an existing index is reused.
        reads: path to a FASTA/FASTQ file of reads to classify.
        output_dir: folder to write all pipeline outputs into (created if
            missing).
        k: k-mer length. Must match the index if index is given.
        index: path to an existing pickled k-mer index. If given, the build
            stage is skipped and this index is reused as-is.
        source: sample source label (e.g. "SRA"), passed through to
            classify_file's database insert.
        top_n: number of species to plot individually before collapsing the
            rest into "Other".
        skip_plot: if True, skip the plotting stage entirely.
        min_quality: if given, drop reads whose mean Phred quality score is
            below this threshold before classifying (FASTQ input only - see
            classify_reads.classify_file / qc.filter_low_quality_reads).
        redistribute: if True, also produce a proportional vote-share
            abundance CSV alongside the winner-take-all one - off by
            default, see classify_reads.classify_file for why.
        db: path to a classification database to add this sample's results
            into, instead of always creating a fresh output_dir/
            classifications.db. Lets multiple pipeline runs (one per reads
            file) accumulate into the same database - pass the same --db
            path on every run in a series, including the first, since
            database.create_database() creates the file if it doesn't
            exist yet. Each run's sample_id (derived from the reads
            filename, see classify_reads.classify_file) must stay distinct
            from any sample already in that database, or its rows get
            replaced (see insert_sample_results).

    Returns:
        dict: paths to the outputs actually produced. Keys: "index",
            "classifications_csv", "database", "diversity_csv",
            "species_abundance_csv", and "plot" (omitted if skip_plot is
            True). species_abundance_csv is get_abundance()'s per-(sample,
            species) breakdown with an added "gc_warning" column (NaN
            unless genome_dir was given and that species was flagged as a
            GC outlier - see qc.gc_outlier_warnings). diversity_csv is
            diversity_by_sample()'s per-sample richness/Shannon output,
            merged with get_classification_totals()'s total/classified/
            unclassified read counts and unclassified_percent.
    """
    if not index and not genome_dir:
        raise ValueError("Either --index or --genome-dir must be provided.")

    os.makedirs(output_dir, exist_ok=True)
    outputs = {}
    gc_warnings = {}
    genomes = None

    if genome_dir:
        click.echo("Checking reference genome GC content...")
        try:
            genomes = load_all_genomes(genome_dir)
        except Exception as e:
            raise RuntimeError(f"Stage 1 (load genomes) failed: {e}") from e
        gc_by_species = {species: gc_content(seq) for species, seq in genomes.items()}
        gc_warnings = gc_outlier_warnings(gc_by_species)
        if gc_warnings:
            for species, message in gc_warnings.items():
                click.echo(f"  Warning ({species}): {message}")
        else:
            click.echo("  No GC outliers detected in reference set.")
    else:
        click.echo("Skipping GC outlier check (no --genome-dir given).")

    if index:
        if not os.path.exists(index):
            raise ValueError(f"--index path does not exist: {index}")
        index_path = index
        click.echo(f"Reusing existing index: {index_path}")
    else:
        click.echo("Stage 1/4: building k-mer index...")
        if genomes is None:
            raise ValueError("--genome-dir is required when --index is not given.")
        try:
            kmer_index = build_kmer_index(genomes, k)
            index_path = os.path.join(output_dir, "kmer_index.pkl")
            save_index(index_path, kmer_index, k)
        except Exception as e:
            raise RuntimeError(f"Stage 1 (build index) failed: {e}") from e
        click.echo(f"Built index with {len(kmer_index)} k-mers, saved to {index_path}")
    outputs["index"] = index_path

    click.echo("Stage 2/4: classifying reads...")
    csv_path = os.path.join(output_dir, "classifications.csv")
    db_path = db if db else os.path.join(output_dir, "classifications.db")
    source_name: str = source if source is not None else ""
    sample_id = derive_sample_id(reads)
    try:
        classify_file(index_path, reads, k, csv_path, db_path=db_path, source=source_name,
                      min_quality=min_quality, redistribute=redistribute)
    except Exception as e:
        raise RuntimeError(f"Stage 2 (classify reads) failed: {e}") from e
    outputs["classifications_csv"] = csv_path
    outputs["database"] = db_path

    click.echo("Stage 3/4: computing diversity metrics...")
    diversity_csv_path = os.path.join(output_dir, "diversity_report.csv")
    species_abundance_path = os.path.join(output_dir, "species_abundance.csv")
    try:
        # get_abundance/get_classification_totals return every sample
        # currently in db_path, not just the one just classified - db_path
        # may be a --db shared with other pipeline runs. Filter down to this
        # run's own sample_id so this run's report/plot describe only this
        # sample, regardless of who else's data lives in the same database.
        abundance_df = get_abundance(db_path)
        abundance_df = abundance_df[abundance_df["sample_id"] == sample_id].reset_index(drop=True)
        abundance_df["gc_warning"] = abundance_df["best_match"].map(gc_warnings)
        abundance_df.to_csv(species_abundance_path, index=False)

        diversity_df = diversity_by_sample(abundance_df)
        totals_df = get_classification_totals(db_path)
        totals_df = totals_df[totals_df["sample_id"] == sample_id].reset_index(drop=True)
        diversity_df = totals_df.merge(diversity_df, on="sample_id", how="left")
        diversity_df[["species_richness", "shannon_diversity"]] = \
            diversity_df[["species_richness", "shannon_diversity"]].fillna(0)
        diversity_df.to_csv(diversity_csv_path, index=False)
    except Exception as e:
        raise RuntimeError(f"Stage 3 (diversity report) failed: {e}") from e
    outputs["diversity_csv"] = diversity_csv_path
    outputs["species_abundance_csv"] = species_abundance_path
    click.echo(diversity_df.to_string(index=False))

    if not skip_plot:
        click.echo("Stage 4/4: plotting species abundance...")
        plot_path = os.path.join(output_dir, "abundance_plot.png")
        try:
            plot_sample_summary(abundance_df, diversity_df.iloc[0], top_n=top_n, output_path=plot_path)
        except Exception as e:
            raise RuntimeError(f"Stage 4 (plot) failed: {e}") from e
        outputs["plot"] = plot_path

    return outputs


@click.command()
@click.option("--genome-dir", default=None,
              help="Folder of reference genomes. Required unless --index is given; also enables the GC-content "
                   "outlier check even when --index is given (skipped if omitted).")
@click.option("--reads", required=True, help="Path to a FASTA/FASTQ file of reads to classify.")
@click.option("--output-dir", required=True, help="Folder to write all pipeline outputs into.")
@click.option("--k", type=int, default=21, help="K-mer length (default: 21).")
@click.option("--index", default=None,
              help="Path to an existing pickled k-mer index. If given, skips the build stage and reuses it.")
@click.option("--source", default=None, help="Sample source label (e.g. 'SRA', 'local') for the database entry.")
@click.option("--top-n", type=int, default=10,
              help="Number of species to plot individually before collapsing the rest into 'Other' (default: 10).")
@click.option("--skip-plot", is_flag=True, default=False, help="Skip the plotting stage.")
@click.option("--min-quality", type=float, default=None,
              help="Minimum mean Phred quality score to keep a read (FASTQ input only). "
                   "Reads below this are dropped before classification.")
@click.option("--redistribute", is_flag=True, default=False,
              help="Also produce a proportional vote-share abundance CSV alongside the "
                   "winner-take-all one. Off by default - not recommended for 16S amplicon "
                   "data, see README 'Validated accuracy'.")
@click.option("--db", default=None,
              help="Path to a classification database to add this sample into, instead of "
                   "creating a fresh output-dir/classifications.db. Use the same --db path "
                   "across multiple runs (different --reads each time) to accumulate samples "
                   "into one database - it's created if it doesn't exist yet. Each run's "
                   "sample_id is derived from the reads filename, so two different reads "
                   "files that produce the same sample_id will overwrite each other in the db.")
def main(genome_dir, reads, output_dir, k, index, source, top_n, skip_plot, min_quality, redistribute, db):
    """Run the full build -> classify -> diversity -> plot pipeline in one command."""
    try:
        outputs = run_pipeline(genome_dir, reads, output_dir, k=k, index=index,
                                source=source, top_n=top_n, skip_plot=skip_plot,
                                min_quality=min_quality, redistribute=redistribute, db=db)
    except (ValueError, RuntimeError) as e:
        raise click.ClickException(str(e))

    click.echo("Pipeline complete. Outputs:")
    for name, path in outputs.items():
        click.echo(f"  {name}: {path}")


if __name__ == "__main__":
    main()
