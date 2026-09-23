#!/bin/env bash
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
