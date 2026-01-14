#!/bin/bash
#SBATCH --job-name=demultiplex_full_v2
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=10:00:00
#SBATCH --error=../slurm_log/%x_%A.err
#SBATCH --output=../slurm_log/%x_%A.out

# Load environment
module purge
module load anaconda3/2024.10
conda activate scseq

cd ../rawdata
mkdir -p r1r3_welldemultiplex_full_v2/r1 r1r3_welldemultiplex_full_v2/r3r1
mkdir -p r1r4_welldemultiplex_full_v2/r1 r1r4_welldemultiplex_full_v2/r4

#  regular 3' adapter, not 3' anchored adapter
srun cutadapt -e 1 \
              -a file:/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/platebarcodes/well_demultiplex_v2.fasta \
              -o ./r1r3_welldemultiplex_full_v2/r1/{name}.1.fastq.gz \
              -p ./r1r3_welldemultiplex_full_v2/r3r1/{name}.31.fastq.gz \
              /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/3380__Merge-sample-PASABC-from-flowcell-HTNFLDRX5-on-2025-05-03_Read_1_passed_filter.fastq.gz \
              /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/trimmed_read3_and_1.fastq.gz &

srun cutadapt -e 1 \
              -a file:/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/platebarcodes/well_demultiplex_v2.fasta \
              -o ./r1r4_welldemultiplex_full_v2/r1/{name}.1.fastq.gz \
              -p ./r1r4_welldemultiplex_full_v2/r4/{name}.4.fastq.gz \
              /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/3380__Merge-sample-PASABC-from-flowcell-HTNFLDRX5-on-2025-05-03_Read_1_passed_filter.fastq.gz \
              /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/3380__Merge-sample-PASABC-from-flowcell-HTNFLDRX5-on-2025-05-03_Read_4_passed_filter.fastq.gz &

wait

