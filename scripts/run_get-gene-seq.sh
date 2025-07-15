#! /bin/bash -login

#SBATCH -D /home/baumlerc/2025-sebastian-paper  # Working directory for the job
#SBATCH -o ./logs/get_gene.%j.out              # Standard output file
#SBATCH -e ./logs/get_gene.%j.err              # Standard error file
#SBATCH -p bmh                             # Partition to submit to
#SBATCH -J get_gene                        # Job name
#SBATCH -t 8:00:00                       # Time limit (8 hours)
#SBATCH -N 1
#SBATCH -n 1                                # Number of tasks
#SBATCH -c 1                                # Number of CPU cores per task
#SBATCH --mem=8G                           # Memory per node
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
    echo "Error: No input csv specified. Please provide one."
    exit 1
fi

file="$1"

sed 1d $file | while IFS=, read -r one two thr four fiv; do 
  echo "Finding sequences for $one of $thr in the $fiv NCBI database"
  python ./get-gene-seq.py --genes $one --organism $thr --extract -db $fiv --all-annotation
done
