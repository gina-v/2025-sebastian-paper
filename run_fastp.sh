#! /bin/bash -login

#SBATCH -D /home/baumlerc/2025-sebastian-paper  # Working directory for the job
#SBATCH -o ./logs/fastp.%j.out              # Standard output file
#SBATCH -e ./logs/fastp.%j.err              # Standard error file
#SBATCH -p high2                            # Partition to submit to
#SBATCH -J fastp                       # Job name
#SBATCH -t 8:00:00                       # Time limit (7 days and 0 hours)
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 1                                # Number of tasks
#SBATCH -c 1                                # Number of CPU cores per task
#SBATCH --mem=64G                           # Memory per node
#SBATCH --mail-type=ALL                     # Send email on all job events
#SBATCH --mail-user=ccbaumler@ucdavis.edu   # Email address for notifications

# Fail on weird errors
set -e
set -x

# Usage: sbatch ./run_fastp.sh /path/to/fastq_dir [optional_output_dir]
set -euo pipefail

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
OUTDIR="${2:-$DIR/fastp_output}"
mkdir -p "$OUTDIR"

# Find all FASTQ(.gz) files and store in array
FILES=($(find "$DIR" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fastq.gz" \)))

# Build associative array of R1 and R2
declare -A R1_FILES
declare -A R2_FILES

# Match common R1/R2 patterns
for f in "${FILES[@]}"; do
    fname=$(basename "$f")
    if [[ "$fname" =~ (_R?1_.*|_1).fastq(.gz)?$ ]]; then
        sample=$(echo "$fname" | sed -E 's/_R?1_|_1.*.fastq(.gz)?//;s/\.fastq(.gz)?$//')
        R1_FILES["$sample"]="$f"
    elif [[ "$fname" =~ (_R?2_.*|_2).fastq(.gz)?$ ]]; then
        sample=$(echo "$fname" | sed -E 's/_R?2_|_2.*.fastq(.gz)?//;s/\.fastq(.gz)?$//')
        R2_FILES["$sample"]="$f"
    fi
done

# Run fastp for matched pairs
for sample in "${!R1_FILES[@]}"; do
    if [[ -n "${R2_FILES[$sample]:-}" ]]; then
        R1="${R1_FILES[$sample]}"
        R2="${R2_FILES[$sample]}"
        base=$(basename "$sample")

        echo "Running fastp on $R1 and $R2..."
        fastp \
            -i "$R1" \
            -I "$R2" \
            -o "$OUTDIR/${base}_R1_fastp.fq.gz" \
            -O "$OUTDIR/${base}_R2_fastp.fq.gz" \
            -h "$OUTDIR/${base}_fastp.html" \
            -j "$OUTDIR/${base}_fastp.json"
    fi
done

