# run all in terminal
#-------------- init
export AF_SAMPLE_DIR=/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/af_workdir
mkdir $AF_SAMPLE_DIR
cd $AF_SAMPLE_DIR

export ALEVIN_FRY_HOME="$AF_SAMPLE_DIR/af_home"
simpleaf set-paths
ulimit -n 2048

#--------------- index

simpleaf index \
--use-piscem \
--output /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/genome/PA14 \
--threads 8 \
--ref-seq /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/genome/PA14/transcripts.fa --overwrite

# simpleaf index \
# --use-piscem \
# --output /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/genome/PA14/spliceu \
# --threads 8 \
# --fasta /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/genome/PA14/GCA_000014625.1_ASM1462v1_genomic.fna \
# --gtf /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/genome/PA14/genomic.gtf \
# --ref-type spliceu

#--------------- chemistry
simpleaf chemistry add -n m3seq -g "1{b[16]u[8]x:}2{r:}" -e rc --version 0.1.0

# looks very complicated, give up for now....
#--------------- quantify

simpleaf quant \
--reads1 /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex/PA14+SA.31.fastq.gz \
--reads2 /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex/PA14+SA.4.fastq.gz \
--threads 16 \
--knee \
--index /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/genome/PA14/index \
--chemistry "1{b[16]u[8]x:}2{r:}" --resolution cr-like \
--anndata-out --use-piscem \
--t2g-map /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/genome/PA14/t2g_simpleaf2.txt \
--output $AF_SAMPLE_DIR/PA14+SA_PA14

simpleaf quant \
--reads1 /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex/PA14+SA.31.fastq.gz \
--reads2 /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex/PA14+SA.4.fastq.gz \
--threads 16 \
--unfiltered-pl /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/whitelist_umitools.txt \
--index /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/genome/PA14/index \
--chemistry "1{b[16]u[8]x:}2{r:}" --resolution cr-like \
--anndata-out --use-piscem \
--t2g-map /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/genome/PA14/t2g_simpleaf2.txt \
--output $AF_SAMPLE_DIR/PA14+SA_PA14_umi_rc --expected-ori rc

simpleaf quant \
--reads1 /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex/PA14+SA.31.fastq.gz \
--reads2 /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex/PA14+SA.4.fastq.gz \
--threads 16 \
--unfiltered-pl /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
--index /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/genome/PA14/index \
--chemistry "1{b[16]u[8]x:}2{r:}" --resolution cr-like \
--anndata-out --use-piscem \
--t2g-map /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/genome/PA14/t2g_simpleaf2.txt \
--output $AF_SAMPLE_DIR/PA14+SA_PA14_737K




simpleaf quant --reads1 ../rawdata/condition_demultiplex_full_v2/PA14+SA.31.fastq.gz --reads2 ../rawdata/condition_demultiplex_full_v2/PA14+SA.4.fastq.gz --threads 8 --unfiltered-pl ../input/whitelist/rc_737K-cratac-v1.txt --index ../input/genome/PA14/alevin/index --chemistry "1{b[16]u[8]x:}2{r:}" --anndata-out --no-piscem --t2g-map ../input/genome/PA14/t2g_simpleaf2.txt --output ../output/alevin/PA14/v2/PA14+SA --resolution cr-like