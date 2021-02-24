#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 12G
#SBATCH --time=96:00:00
#SBATCH --partition=mrcieu
set -euo pipefail

module load languages/r/3.6.0
module load apps/qctool/2.0rc4
module load apps/bgen/1.1.6

# run simulation
Rscript sim.R
