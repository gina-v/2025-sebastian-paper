#! /usr/bin/env python

import polars as pl
import argparse
import os

def main():
    p = argparse.ArgumentParser(description="Merge mapping and coverage stats by contig.")
    p.add_argument(
        "run_accession",
        help="Run accession used in filename (e.g. SRR1211428) without suffix. Will look for *_mapping_stats.tsv and *_average_coverage_stats.tsv"
    )
    p.add_argument(
        "protein_accession",
        help="Protein accession used in filename (e.g. NP_462265.1) without suffix. Will look for *_mapping_stats.tsv and *_average_coverage_stats.tsv"
    )
    p.add_argument(
        "-d", "--directory", default=None,
        help="The directory to look for *_mapping_stats.tsv and *_average_coverage_stats.tsv"
    )
    p.add_argument(
        "-o", "--output", default=None,
        help="Output file path. If not provided, prints to stdout."
    )

    args = p.parse_args()
    run_a = args.run_accession
    pro_a = args.protein_accession

    if args.directory:
        try:
            base = os.path.join(args.directory, run_a + ".x." + pro_a)
        except ValueError:
            raise ValueError("Basename must be of the form SRRxxxxxxx.x.NP_xxxxxx.x")
    else:
        print("No directory argument used. Search current directory...")
        try:
            base = run_a + ".x." + pro_a
        except ValueError:
            raise ValueError("Basename must be of the form SRRxxxxxxx.x.NP_xxxxxx.x")
 
    map_file = f"{base}_mapping_stats.tsv"
    cov_file = f"{base}_average_coverage_stats.tsv"

    if not os.path.isfile(map_file):
        raise FileNotFoundError(f"Mapping stats file not found: {map_file}")
    if not os.path.isfile(cov_file):
        raise FileNotFoundError(f"Coverage stats file not found: {cov_file}")

    mapping_stats = pl.read_csv(map_file, separator="\t")
    coverage_stats = pl.read_csv(cov_file, separator="\t")

    df = mapping_stats.join(coverage_stats, on="contig", how="inner")

    df = df.with_columns([
        pl.lit(run_a).alias("run"),
        pl.lit(pro_a).alias("gene")
    ])

    df = df.select(["contig", "run", "gene"] + [col for col in df.columns if col not in {"contig", "run", "gene"}])

    if df.is_empty():
        # Create a new DataFrame with one row
        new_row = pl.DataFrame({
            "contig": ["None"],
            "run": [run_a],
            "gene": [pro_a],
            **{col: [0] for col in df.columns if col not in {"contig", "run", "gene"}}
        })

        # Assign the new DataFrame
        df = new_row


    if args.output:
        df.write_csv(args.output, separator=",")
        print(f"Output written to {args.output}")
    else:
        print(df)


if __name__ == "__main__":
    main()

