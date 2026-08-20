import argparse
from database import get_abundance, get_classification_totals
from diversity import diversity_by_sample
from visualization import plot_multi_sample_abundance


def main():
    parser = argparse.ArgumentParser(description="Report per-sample species richness and Shannon diversity from a classification database.")
    parser.add_argument("--db", required=True, help="Path to a SQLite database populated by classify_reads.py (via --db)")
    parser.add_argument("--output", required=False, help="Path to save the diversity report (.csv)")
    parser.add_argument("--plot", required=False, help="Path to save a species-abundance-by-sample chart (.png)")
    parser.add_argument("--top-n", type=int, default=10, help="Number of species to plot individually before collapsing the rest into 'Other' (default: 10)")

    args = parser.parse_args()

    # Queried once and reused: the plot needs the same per-(sample, species)
    # breakdown the diversity metrics are derived from, and re-querying it
    # means a second full aggregation over every classified read.
    abundance_df = get_abundance(args.db)
    report_df = diversity_by_sample(abundance_df)
    print(report_df.to_string(index=False))

    if args.output:
        report_df.to_csv(args.output, index=False)
        print(f"Diversity report saved to {args.output}")

    if args.plot:
        totals_df = get_classification_totals(args.db)
        unclassified_percent = totals_df.set_index("sample_id")["unclassified_percent"]
        plot_multi_sample_abundance(abundance_df, unclassified_percent, top_n=args.top_n, output_path=args.plot)
        print(f"Abundance plot saved to {args.plot}")


if __name__ == "__main__":
    main()
