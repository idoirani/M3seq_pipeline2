#!/bin/bash
#SBATCH --job-name=combash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=2:30:00
#SBATCH --error=../slurm_log/%x_%A.err
#SBATCH --output=../slurm_log/%x_%A.out

# Load environment
module purge
module load anaconda3/2024.10
conda activate scseq

export PIGZ_OPTS="-p $SLURM_CPUS_PER_TASK"

# Input: trimmed FASTQ.GZ files from previous cutadapt run
TRIM1=../rawdata/trimmed_read1.fastq.gz
TRIM3=../rawdata/trimmed_read3.fastq.gz
OUT=../rawdata/trimmed_read3_and_1.fastq.gz

# Step: Combine trimmed reads
python combine_fast.py "$TRIM3" "$TRIM1" "$OUT"
