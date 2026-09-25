#!/bin/env bash
# =============================================================================
#  Slurm array job for the L50 run set: one task per row of DF.csv.
#
#  Submit from the simulation directory so that $SLURM_SUBMIT_DIR points at it:
#      sbatch --chdir=2D ../slurm/submit.sh
#
#  Each task gets one GPU and runs 2D.jl with SLURM_ARRAY_TASK_ID as the
#  simulation index; set DATA_DIR in your environment (or edit
#  InputParameters.jl) to choose where snapshots are written.
#
#  Reconstructed from the job generator used for the published runs. The
#  partition and constraint name one particular cluster — change them for
#  wherever you are running.
#
#    --array=1-13%20   13 simulations, at most 20 running at once
#    --time            12 h, which is what t_fin = 150000 was sized against
#    --mem=3000        host memory; the fields live on the GPU
#    --constraint      A100 40 GB or 80 GB (the 5008² Float64 grid needs it)
# =============================================================================
#SBATCH --array=1-13%20
#SBATCH --partition=private-kruse-gpu,shared-gpu
#SBATCH --time=0-12:00:00
#SBATCH --output=%J.out
#SBATCH --mem=3000
#SBATCH --gpus=1
#SBATCH --constraint=nvidia_a100-pcie-40gb|nvidia_a100_80gb_pcie

module load Julia

cd "$SLURM_SUBMIT_DIR"
srun julia --optimize=3 2D.jl
