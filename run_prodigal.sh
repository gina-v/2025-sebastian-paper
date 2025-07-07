#! /bin/bash -login

#SBATCH -D /home/baumlerc/2025-sebastian-paper  # Working directory for the job
#SBATCH -o ./logs/prodigal.%j.out              # Standard output file
#SBATCH -e ./logs/prodigal.%j.err              # Standard error file
#SBATCH -p bmh                             # Partition to submit to
#SBATCH -J gene_annotation                       # Job name
#SBATCH -t 8:00:00                       # Time limit (8 hours)
#SBATCH -N 1                                # Number of nodes
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G                           # Memory per node
#SBATCH --mail-type=ALL                     # Send email on all job events
#SBATCH --mail-user=ccbaumler@ucdavis.edu   # Email address for notifications

# Usage: sbatch ./run_prodigal.sh /path/to/metahit_dir [optional_output_dir]
set -euo pipefail
set -x

# Initialize conda if needed
 conda_base=$(conda info --base)
 . ${conda_base}/etc/profile.d/conda.sh

# Activate relevant environment if needed
conda activate sebastian

# Check if any directory arguments are provided; exit if none are found
if [ "$#" -eq 0 ]; then
    echo "Error: No directories specified. Please provide one or more directories to check."
    exit 1
fi

DIR="${1:-.}"
OUTDIR="${2:-$DIR/prodigal_output}"
mkdir -p "$OUTDIR"

declare -A CONTIGS

while IFS= read -r -d '' f; do
    sample=$(dirname "$f")
    CONTIGS["$sample"]="$f"
done < <(find "$DIR" -maxdepth 2 -type f -name "final.contigs.fa" -print0)

for dir in "${!CONTIGS[@]}"; do
    echo "Directory: $sample"
    echo "Contig file: ${CONTIGS[$sample]}"
done

# https://genomicsaotearoa.github.io/hts_workshop_mpi/level2/42_annotation_prodigal/#predicting-protein-coding-regions
for sample in "${!CONTIGS[@]}"; do
    echo "Running prodigal on $sample..."
    base=$(basename $sample)
    echo "for $base"
    echo "${CONTIGS[$sample]}"

    # Metagenome mode, # Treat runs of N as masked sequence; dont build genes across them.
    prodigal \
        -p meta \
        -m \
        -i "${CONTIGS[$sample]}" \
        -d "$OUTDIR${base}.prod.fna" \
        -a "$OUTDIR${base}.prod.faa" \
        -o "$OUTDIR${base}.prod.gbk"
done
