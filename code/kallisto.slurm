## have to run in terminal with della9
cd /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/genome/PA14
mkdir kallisto
cd kallisto

kb ref -i index.idx -g t2g.txt -f1 transcriptome.fa ../GCA_000014625.1_ASM1462v1_genomic.fna ../genomic.gtf

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 -o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/PA14+SA/PA14/kallisto/rc --overwrite -w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex/PA14+SA.31.fastq.gz /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex/PA14+SA.4.fastq.gz --h5ad --report 

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 -o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/PA14+SA/PA14/kallisto/umitools --overwrite -w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/whitelist_umitools.txt /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex/PA14+SA.31.fastq.gz /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex/PA14+SA.4.fastq.gz --h5ad --report

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 -o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/PA14+BC/PA14/kallisto/rc --overwrite -w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex/PA14+BC.31.fastq.gz /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex/PA14+BC.4.fastq.gz --h5ad --report

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 -o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/PA14+SA/PA14/kallisto/737K --overwrite -w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex/PA14+SA.31.fastq.gz /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex/PA14+SA.4.fastq.gz --h5ad --report --strand unstranded

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 -o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/PA14+SA/PA14/kallisto/fastp --overwrite -w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA.31.fastp2.fastq.gz /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA.4.fastp2.fastq.gz --h5ad --report --strand unstranded

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 -o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/PA14+SA/PA14/kallisto/old_PA --overwrite -w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt /scratch/gpfs/ziyangc/PASA/ScRNAseq/old/data/rawdata/PA.trimmed.fastq.gz /scratch/gpfs/ziyangc/PASA/ScRNAseq/old/data/rawdata/PA.read4.fastq.gz --h5ad --report --strand unstranded

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 -o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/PA14+SA/PA14/kallisto/v2_PA --overwrite -w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PA14+SA.31.fastq.gz /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PA14+SA.4.fastq.gz --h5ad --report --strand unstranded

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 -o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/PA14+SA/PA14/kallisto/v2_PA --overwrite -w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PAO1+SA.31.fastq.gz /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PAO1+SA.4.fastq.gz --h5ad --report --strand unstranded

#--------------------------------------

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/test_combination\PA14_index\v2\PA14+SA --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PA14+SA.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PA14+SA.4.fastq.gz \
--h5ad --report --strand unstranded

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/test_combination\PA14_index\v2\PA14+BC --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PA14+BC.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PA14+BC.4.fastq.gz \
--h5ad --report --strand unstranded

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/test_combination\PA14_index\v2\PAO1+BC --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PAO1+BC.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PAO1+BC.4.fastq.gz \
--h5ad --report --strand unstranded

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/kallisto_combination/PA14_index/v2/PAO1+SA --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PAO1+SA.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PAO1+SA.4.fastq.gz \
--h5ad --report --strand unstranded


kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/kallisto_combination/PA14_index/v1/PAO1+SA --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PAO1+SA.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PAO1+SA.4.fastq.gz \
--h5ad --report --strand unstranded

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/kallisto_combination/PA14_index/v1/PAO1+BC --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PAO1+BC.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PAO1+BC.4.fastq.gz \
--h5ad --report --strand unstranded

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/kallisto_combination/PA14_index/v1/PA14+BC --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+BC.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+BC.4.fastq.gz \
--h5ad --report --strand unstranded

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/kallisto_combination/PA14_index/v1/PA14+SA --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA.4.fastq.gz \
--h5ad --report --strand unstranded

#---------------------------------- PAO1

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/kallisto_combination/PAO1_index/v1/PA14+SA --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA.4.fastq.gz \
--h5ad --report --strand unstranded

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/kallisto_combination/PAO1_index/v1/PA14+BC --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+BC.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+BC.4.fastq.gz \
--h5ad --report --strand unstranded

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/kallisto_combination/PAO1_index/v1/PAO1+BC --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PAO1+BC.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PAO1+BC.4.fastq.gz \
--h5ad --report --strand unstranded

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/kallisto_combination/PAO1_index/v1/PAO1+SA --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PAO1+SA.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PAO1+SA.4.fastq.gz \
--h5ad --report --strand unstranded

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/kallisto_combination/PAO1_index/v2/PA14+SA --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PA14+SA.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PA14+SA.4.fastq.gz \
--h5ad --report --strand unstranded

#--------------------------------- USA300

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/kallisto_combination/USA300_index/v1/PA14+SA --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA.4.fastq.gz \
--h5ad --report --strand unstranded

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/kallisto_combination/USA300_index/v2/PA14+SA --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PA14+SA.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PA14+SA.4.fastq.gz \
--h5ad --report --strand unstranded


#------------------ PA14-seq

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/kallisto_combination/PA14_seq/v1/PA14+SA --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA.4.fastq.gz \
--h5ad --report --strand unstranded


kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/kallisto_combination/PA14_seq/v2/PA14+SA --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PA14+SA.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PA14+SA.4.fastq.gz \
--h5ad --report --strand unstranded

#------------------- BC

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/kallisto_combination/BC_index/v1/PA14+SA --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA.4.fastq.gz \
--h5ad --report --strand unstranded

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/kallisto_combination/BC_index/v1/PA14+BC --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+BC.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+BC.4.fastq.gz \
--h5ad --report --strand unstranded

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/kallisto_combination/BC_index/v2/PA14+SA --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PA14+SA.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PA14+SA.4.fastq.gz \
--h5ad --report --strand unstranded

kb count -i index.idx -g t2g.txt -x 0,0,16:0,16,24:1,0,0 \
-o /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/output/kallisto_combination/BC_index/v2/PA14+BC --overwrite \
-w /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/input/whitelist/rc_737K-cratac-v1.txt \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PA14+BC.31.fastq.gz \
/scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full_v2/PA14+BC.4.fastq.gz \
--h5ad --report --strand unstranded