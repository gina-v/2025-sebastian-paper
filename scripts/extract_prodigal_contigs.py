#! /usr/bin/env python

import os
import argparse
import csv
import sys
from collections import defaultdict

def main():
    p = argparse.ArgumentParser(description="Group the seq_filename and sseq_id for the parsed BLAST output and map the reads to the contigs")

    p.add_argument('input', help='The parsed BLAST csv file')
    p.add_argument('-o', '--output', help='The output file for mapping reads')
    p.add_argument('-v', '--verbose', action='store_true')
    p.add_argument('-d', '--prodigal-dir', help='The file path to the prodigal directory')
    p.add_argument('-f', '--prodigal-filetype', help='The desired contig file type from prodigal')
    p.add_argument('--outdir', help='The path to the contig file output directory')
    p.add_argument('--possible-samples', nargs='+', default=[],
                   help='List of all possible sample seq_filenames to create empty files if missing')

    args = p.parse_args()

    with open(args.input, newline='') as fp:
        reader = csv.DictReader(fp)

        seq_dict = defaultdict(list)
        first_qseqid = None

        for row in reader:
            qseqid = row.get('qseqid')
            if first_qseqid is None:
                first_qseqid = qseqid
            else:
                assert qseqid == first_qseqid, f"Error: Mismatch in 'qseqid' column.\nExpected {first_qseqid}, got {qseqid}"

            seq_filename = row.get('seq_filename', '')
            sseq_id = row.get('sseq_id', '')
            assert seq_filename, "Error: Empty 'seq_filename' found"
            seq_dict[seq_filename].append(sseq_id)

    if args.verbose:
        print(f"\nAll rows matched qseqid: {first_qseqid}")
        print("\nGrouped sseq_id by seq_filename:\n")
        for k, v in seq_dict.items():
            print(f"{k}: {', '.join(v)}")

    for seq_file, contigs in seq_dict.items():
        found = 0
        total = len(contigs)
        write_record = False

        full_prod_path = os.path.join(args.prodigal_dir, seq_file.strip().split(".")[0] + ".prod." + args.prodigal_filetype)
        full_contig_path = os.path.join(args.outdir, seq_file.strip().split(".")[0] + ".x." + qseqid + ".contiglist." + args.prodigal_filetype)
        if args.verbose: print(f"\nRetrieving {full_prod_path}")

        with open(full_prod_path) as fin, open(full_contig_path, 'w') as fout:
            for line in fin:
                if line.startswith('>'):
                    contig_id = line[1:].strip().split()[0]
                    if contig_id in contigs:
                        write_record = True
                        found += 1
                        fout.write(line)
                    else:
                        write_record = False
                elif write_record:
                    fout.write(line)

        print(f"Finished. Found {found} of {total} contigs for {seq_file}.\nFor the BLAST'd sequence {qseqid}:\nSample {seq_file}: {', '.join(contigs)}")

        if args.verbose: print(f"Writing to {full_contig_path}\n")

    expected_samples_set = set(args.possible_samples)
    used_samples_set = set(seq_dict.keys())
    missing_samples = expected_samples_set - used_samples_set

    for missing_sample in missing_samples:
        empty_contig_path = os.path.join(args.outdir, missing_sample.strip().split(".")[0] + ".x." + first_qseqid + ".contiglist." + args.prodigal_filetype)
        if not os.path.exists(empty_contig_path):
            if args.verbose:
                print(f"Creating empty contiglist file for missing sample: {empty_contig_path}")
            with open(empty_contig_path, "w") as fout:
                pass  # create empty file

    if args.output:
        with open(args.output, 'w', newline='') as fout:
            writer = csv.writer(fout)
            writer.writerow(['seq_filename', 'sseq_ids'])
            for k, v in seq_dict.items():
                writer.writerow([k, ','.join(v)])

if __name__ == "__main__":
    sys.exit(main())
