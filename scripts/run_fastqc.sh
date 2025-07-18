#! /bin/bash -login

#SBATCH -D /home/baumlerc/2025-sebastian-paper  # Working directory for the job
#SBATCH -o ./logs/fastqc.%j.out              # Standard output file
#SBATCH -e ./logs/fastqc.%j.err              # Standard error file
#SBATCH -p high2                            # Partition to submit to
#SBATCH -J fastqc                       # Job name
#SBATCH -t 8:00:00                       # Time limit (8 hours)
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 1                                # Number of tasks
#SBATCH -c 1                                # Number of CPU cores per task
#SBATCH --mem=64G                           # Memory per node
#SBATCH --mail-type=ALL                     # Send email on all job events
#SBATCH --mail-user=ccbaumler@ucdavis.edu   # Email address for notifications

# Usage: ./run_fastp.sh /path/to/fastq_dir [optional_output_dir]
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
OUTDIR="${2:-$INDIR/fastqc_output}"
mkdir -p "$OUTDIR"

# Find all FASTP.FQ(.gz) files and store in array
FILES=($(find "$DIR" -maxdepth 1 -type f \( -name "*fastp.fq" -o -name "*fastp.fq.gz" \)))

# Build associative array of R1 and R2
declare -A R1_FILES
declare -A R2_FILES

# Match common R1/R2 patterns
for f in "${FILES[@]}"; do
    fname=$(basename "$f")

    if [[ "$fname" == *_R1_fastp.fq.gz || "$fname" == *_1_fastp.fq.gz ]]; then
        sample=$(echo "$fname" | sed -E 's/_R?1_fastp\.fq\.gz$//; s/_1_fastp\.fq\.gz$//')
        R1_FILES["$sample"]="$f"

    elif [[ "$fname" == *_R2_fastp.fq.gz || "$fname" == *_2_fastp.fq.gz ]]; then
        sample=$(echo "$fname" | sed -E 's/_R?2_fastp\.fq\.gz$//; s/_2_fastp\.fq\.gz$//')
        R2_FILES["$sample"]="$f"
    fi
done

# Run fastqc for matched pairs
for sample in "${!R1_FILES[@]}"; do
    echo $sample
    if [[ -n "${R2_FILES[$sample]:-}" ]]; then
        R1="${R1_FILES[$sample]}"
        R2="${R2_FILES[$sample]}"
        base=$(basename "$sample")

        mkdir -p $OUTDIR${base}

        echo "Running fastqc on $R1 and $R2..."
        fastqc \
            "$R1" \
            "$R2" \
            --outdir "$OUTDIR${base}"
    fi
done

