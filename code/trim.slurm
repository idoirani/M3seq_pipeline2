#!/bin/bash
#SBATCH --array=1
#SBATCH --job-name=trim         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=3               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=8G         # memory per cpu-core (4G is default)
#SBATCH --time=5:01:00          # total run time limit (HH:MM:SS)
#SBATCH --error=../slurm_log/%x_%A.err
#SBATCH --output=../slurm_log/%x_%A.out

module purge
module load anaconda3/2024.10
conda activate scseq

cd ../rawdata

# Input files
READ1=/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/3380__Merge-sample-PASABC-from-flowcell-HTNFLDRX5-on-2025-05-03_Read_1_passed_filter.fastq.gz
READ3=/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/3380__Merge-sample-PASABC-from-flowcell-HTNFLDRX5-on-2025-05-03_Read_3_Index_Read_passed_filter.fastq.gz

# Temporary trimmed outputs
TRIM1=trimmed_read1.fastq.gz
TRIM3=trimmed_read3.fastq.gz
INDEX=trimmed_read1_platebarcode.fastq.gz

# # Step 1: Trim 3' 18 nt from read1, and 5' 14 nt from read3, also the 5' 8 nt from read1 so get anchored platebarcode for faster demultiplex
cutadapt -u -18 -o "$TRIM1" "$READ1" &
cutadapt -u 14 -o "$TRIM3" "$READ3" &
cutadapt -u 8 -o "$INDEX" "$READ1" &
wait