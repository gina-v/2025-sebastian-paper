#!/usr/bin/env python

import argparse
import csv
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Filter CSV by e-value and split sseqid.")
    parser.add_argument('--cutoff', type=float, help='E-value cutoff. Only rows with e-value < this will be kept.')
    parser.add_argument('input', type=argparse.FileType('r'), help='Input CSV file (or stdin).')
    parser.add_argument('-o', '--output', required=True, help='Output CSV file.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print verbose output to stdout.')
    return parser.parse_args()

def split_sseqid(sseqid):
    if '.faa' in sseqid and '_' in sseqid:
        split_index = sseqid.find('_')
        return sseqid[:split_index], sseqid[split_index+1:]
    else:
        return '', sseqid

def main():
    args = parse_args()

    # Check for header
    sample = args.input.read(2048)
    args.input.seek(0)
    has_header = csv.Sniffer().has_header(sample)
    if not has_header:
        sys.exit(f"Error: Input file does not appear to have a header row.\n{sample}")

    expected_headers = ['qseqid', 'sseqid', 'pident', 'length', 'e_value', 'bit_score']
    reader = csv.reader(sample.splitlines())
    header_row = next(reader)

    missing = [h for h in expected_headers if h not in header_row]
    if missing:
        sys.exit(f"Error: Missing expected columns: {', '.join(missing)}.\nFound: {header_row}\nUse the `-outfmt '10 qseqid sseqid pident length evalue bitscore'` in BLAST CLI \nAnd then on the output use `cat <(printf 'qseqid,sseqid,pident,length,e_value,bit_score') blast_output > blast_output_header`")

    reader = csv.DictReader(args.input)

    # Ensure fieldnames match output dict
    fieldnames = ['qseqid', 'seq_filename', 'sseq_id', 'pident', 'length', 'e_value', 'bit_score']
    if args.verbose:
        print(" | ".join(f"{f}" for f in fieldnames))

    with open(args.output, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            try:
                e_value = float(row['e_value'])

                if args.cutoff is None or e_value < args.cutoff: #The closer to 0 the e-value the better

                    seq_filename, sseq_id = split_sseqid(row['sseqid'])

                    output_row = {
                        'qseqid': row['qseqid'],
                        'seq_filename': seq_filename,
                        'sseq_id': sseq_id,
                        'pident': row['pident'],
                        'length': row['length'],
                        'e_value': row['e_value'],
                        'bit_score': row['bit_score']
                    }

                    if args.verbose:
                        print(" | ".join(f"{k}: {v}" for k, v in output_row.items()))

                    writer.writerow(output_row)

            except ValueError:
                continue  # Skip rows with invalid e_value

if __name__ == "__main__":
    main()
