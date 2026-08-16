import argparse
from database import get_abundance
from diversity import diversity_by_sample


def generate_report(db_path: str):
    """
    Compute per-sample diversity metrics from a classification database.

    Parameters:
        db_path (str): Path to a SQLite database populated by
            classify_reads.py (via --db).

    Returns:
        pd.DataFrame: diversity_by_sample()'s output - one row per sample_id
            with species_richness and shannon_diversity.
    """
    abundance_df = get_abundance(db_path)
    return diversity_by_sample(abundance_df)


def main():
    parser = argparse.ArgumentParser(description="Report per-sample species richness and Shannon diversity from a classification database.")
    parser.add_argument("--db", required=True, help="Path to a SQLite database populated by classify_reads.py (via --db)")
    parser.add_argument("--output", required=False, help="Path to save the diversity report (.csv)")

    args = parser.parse_args()

    report_df = generate_report(args.db)
    print(report_df.to_string(index=False))

    if args.output:
        report_df.to_csv(args.output, index=False)
        print(f"Diversity report saved to {args.output}")


if __name__ == "__main__":
    main()
