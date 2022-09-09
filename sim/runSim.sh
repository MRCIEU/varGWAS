#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem 10700
#SBATCH --time=72:00:00
#SBATCH --partition=cpu,mrcieu
set -euo pipefail

module load apps/singularity/3.8.3

mkdir -p data
singularity exec \
--no-mount home \
--bind /user/home/ml18692/projects:/user/home/ml18692/projects \
--pwd `pwd` \
/user/home/ml18692/projects/varGWAS/sim/vargwas-sim.sif \
Rscript "$@"
