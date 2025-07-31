#!/bin/bash
set -euo pipefail

if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <output_dir> <threads> <sample_metadata_file> <gene_metadata_file>"
  exit 1
fi

o="$1"
threads="$2"
sample_metadata="$3"
gene_metadata="$4"

mkdir -p "${o}/mapped_reads/"

sample_delim="\t"
symbol_delim=","
max_rows=0  # Set to 0 or leave empty to read all rows

# Read SAMPLES array
SAMPLES=()
if [[ "$max_rows" -gt 0 ]]; then
    SAMPLES=($(awk -F"$sample_delim" -v max="$max_rows" '
        NR==1 {for (i=1; i<=NF; i++) if ($i=="run") col=i; next}
        NR>1 && NR<=max+1 {print $col}
    ' "$sample_metadata"))
else
    SAMPLES=($(awk -F"$sample_delim" '
        NR==1 {for (i=1; i<=NF; i++) if ($i=="run") col=i; next}
        {print $col}
    ' "$sample_metadata"))
fi

# Read SYMBOLS array
SYMBOLS=()
if [[ "$max_rows" -gt 0 ]]; then
    SYMBOLS=($(awk -F"$symbol_delim" -v max="$max_rows" '
        NR==1 {for (i=1; i<=NF; i++) if ($i=="symbol") col=i; next}
        NR>1 && NR<=max+1 {print $col}
    ' "$gene_metadata"))
else
    SYMBOLS=($(awk -F"$symbol_delim" '
        NR==1 {for (i=1; i<=NF; i++) if ($i=="symbol") col=i; next}
        {print $col}
    ' "$gene_metadata"))
fi

# Debug output (optional)
echo "Loaded ${#SAMPLES[@]} samples"
echo "Loaded ${#SYMBOLS[@]} symbols"

for sample in "${SAMPLES[@]}"; do
  for symbol in "${SYMBOLS[@]}"; do
    contiglist_csv="${o}/contiglists/${symbol}_contiglist.csv"
    fq1="${o}/fastq/${sample}_1.fastq.gz"
    fq2="${o}/fastq/${sample}_2.fastq.gz"
    contig_list_fna="${o}/contiglists/${sample}.x.${symbol}.contiglist.fna"
    sorted_bam="${o}/mapped_reads/${sample}.x.${symbol}.sorted.bam"

    echo "Running minimap2 for sample=$sample symbol=$symbol"
    minimap2 -t "${threads}" \
      "${contig_list_fna}" \
      "${fq1}" "${fq2}" \
    | samtools sort -@ "${threads}" -o "${sorted_bam}"

    echo "Finished: ${sorted_bam}"
  done
done
