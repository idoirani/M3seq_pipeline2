#!/bin/bash
#SBATCH --job-name=whitelist_distance
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=06:30:00
#SBATCH --error=../slurm_log/%x_%A.err
#SBATCH --output=../slurm_log/%x_%A.out

# Load environment
module purge
module load anaconda3/2024.10
conda activate scseq

bash whitelist.sh
