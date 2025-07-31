#! /usr/bin/env python

import polars as pl
import argparse
import os


def check_file(file):
    if not os.path.isfile(file):
        raise FileNotFoundError(f"Mapping stats file not found: {file}")
    return file

def fix_numeric_col_type(df):
    df = df.with_columns([
        pl.col(col).cast(pl.Float64) 
        for col, dtype in df.schema.items() 
        if dtype.is_numeric()
    ])
    return df

def main():
    p = argparse.ArgumentParser(description="Merge mapping and coverage stats by contig.")
    p.add_argument(
        "-m", "--metadata-file",
        help="The metadata file containing the diagnosis and read_count metrics"
    )
    p.add_argument(
        "-i", "--input-directory", default=None,
        help="The directory that contains the csv files of the merged_stats.py"
    )
    p.add_argument(
        "-o", "--output", default=None,
        help="Output file path. If not provided, prints to stdout."
    )

    args = p.parse_args()
    indir = args.input_directory

    mf = pl.read_csv(args.metadata_file, separator='\t')

    csv_files = [os.path.join(indir, f) 
                 for f in os.listdir(indir) 
                 if f.endswith(".csv")]

    if not csv_files:
        print("No CSV files found in the directory.")
        return

    contigs_df = pl.concat([fix_numeric_col_type(pl.read_csv(check_file(fp))) for fp in csv_files])
    print(contigs_df)

    print(mf)
    mf_less = mf.select(["run", "diagnosis", "read_count"])
    joined_df = contigs_df.join(mf_less, on="run", how="left")
    df = joined_df.select(["contig", "run", "gene", "diagnosis", "read_count"] + [col for col in joined_df.columns if col not in {"contig", "run", "gene", "diagnosis", "read_count"}])

    if args.output:
        df.write_csv(args.output, separator=",")
        print(f"Output written to {args.output}")
    else:
        print(df)

if __name__ == "__main__":
    main()

