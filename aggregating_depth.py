#! /usr/bin/env python

import argparse
from collections import defaultdict

def main():
    p = argparse.ArgumentParser(description="Collapse rows by contig and compute max position and average count.")
    p.add_argument("input", help="Input file with three tab-separated columns: contig, position, count.")
    p.add_argument("output", help="Output file to write the collapsed results.")
    args = p.parse_args()

    positions = defaultdict(list)
    counts = defaultdict(list)

    with open(args.input) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 3:
                continue
            contig, pos, count = parts
            positions[contig].append(int(pos))
            counts[contig].append(int(count))

    with open(args.output, 'w') as out:
        out.write("contig\tlength\taverage_depth\n")
        for contig in sorted(positions):
            max_pos = max(positions[contig])
            avg_count = sum(counts[contig]) / len(counts[contig])
            out.write(f"{contig}\t{max_pos}\t{avg_count:.2f}\n")

if __name__ == "__main__":
    main()
