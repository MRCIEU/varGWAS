#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem 28571
#SBATCH --time=12:00:00
#SBATCH --partition=mrcieu,cpu
set -euo pipefail

module load languages/r/3.6.0
module load apps/qctool/2.0rc4
module load apps/bgen/1.1.6
module load languages/gcc/9.3.0

# run
mkdir -p data
PATH="$PATH":/mnt/storage/home/ml18692/projects/varGWAS/build/bin
PATH="$PATH":/mnt/storage/home/ml18692/apps/osca
Rscript "$@"
