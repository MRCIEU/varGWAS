#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem 28571
#SBATCH --time=12:00:00
#SBATCH --partition=mrcieu,cpu
set -euo pipefail

module load apps/singularity/3.8.3

# run
mkdir -p data
singularity exec \
--no-mount home \
--bind /user/home/ml18692/projects:/user/home/ml18692/projects \
--pwd `pwd` \
/user/home/ml18692/projects/varGWAS/sim/vargwas-sim.sif \
Rscript "$@"