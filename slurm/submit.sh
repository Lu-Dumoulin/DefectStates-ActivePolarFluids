#!/bin/env bash
# =============================================================================
#  Slurm array job: one task per row of a parameter table.
#
#  Submit from the repository root:
#
#      export DF_FILE=../params/DF_9.csv       # which table (default: 2D/DF.csv)
#      export DATA_DIR=/scratch/me/fig9/       # where snapshots go
#      sbatch --chdir=2D --array=1-1680%20 slurm/submit.sh
#
#  --array must match the number of rows in the table. Slurm exports your
#  environment to the job, so DF_FILE and DATA_DIR reach every task, and each
#  task runs the row given by its SLURM_ARRAY_TASK_ID.
#
#  The partition and constraint below name one particular cluster: change them
#  for yours. The time limit suits the L=50 runs; the 2e6 lattice runs
#  (figures 8, 12, 13) need several days.
# =============================================================================
#SBATCH --array=1-13%20
#SBATCH --partition=private-kruse-gpu,shared-gpu
#SBATCH --time=0-12:00:00
#SBATCH --output=%J.out
#SBATCH --mem=3000
#SBATCH --gpus=1
#SBATCH --constraint=nvidia_a100-pcie-40gb|nvidia_a100_80gb_pcie

module load Julia

# The job starts in the directory given by --chdir, i.e. 2D/.
srun julia --project=.. --optimize=3 2D.jl
