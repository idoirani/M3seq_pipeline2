salmon quant -l A \
-i /scratch/gpfs/ziyangc/PASA/BulkRNASeq3/input/genbank/ASM1462/ncbi_dataset/data/GCA_000014625.1/salmon_index \
-r ../rawdata/condition_demultiplex_full/PA14+SA.4.fastq.gz \
-o ../output/salmon/PA14/v1/PA14+SA --writeUnmappedNames


awk '{print $1 " 4:N:0:0"}' unmapped_names.txt > modified_names.txt

seqkit grep -n -f unmapped_names.txt /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA.4.fastq.gz \
-o unmapped_reads_v1.fastq


#----------------

salmon quant -l A \
-i /scratch/gpfs/ziyangc/PASA/BulkRNASeq3/input/genbank/ASM1462/ncbi_dataset/data/GCA_000014625.1/salmon_index \
-r /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA.4.fastp2.fastq.gz \
-o PA14+SA_salmon/PA14 --writeUnmappedNames

cd PA14+SA_salmon/PA14/aux_info/

awk '{print $1 " 4:N:0:0"}' unmapped_names.txt > modified_names.txt

seqkit grep -n -f modified_names.txt /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA.4.fastp2.fastq.gz \
-o unmapped_reads.fastq

gzip unmapped_reads.fastq

cd ../../../

salmon quant -l A \
-i /scratch/gpfs/ziyangc/PASA/BulkRNASeq3/input/genbank/ASM676+25608/bakta/salmon_index \
-r /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA_salmon/PA14/aux_info/unmapped_reads.fastq.gz \
-o PA14+SA_salmon/PA14+PAO1+BC --writeUnmappedNames

cd PA14+SA_salmon/PA14+PAO1+BC/aux_info/

awk '{print $1 " 4:N:0:0"}' unmapped_names.txt > modified_names.txt

seqkit grep -n -f modified_names.txt /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA_salmon/PA14/aux_info/unmapped_reads.fastq.gz \
-o unmapped_reads.fastq

gzip unmapped_reads.fastq

cd ../../../

salmon quant -l A \
-i /scratch/gpfs/ziyangc/PASA/BulkRNASeq2/input/genbank/ASM1346/ncbi_dataset/data/GCA_000013465.1/salmon_index \
-r /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA_salmon/PA14+PAO1+BC/aux_info/unmapped_reads.fastq.gz \
-o PA14+SA_salmon/PA14+PAO1+BC+SA --writeUnmappedNames

cd PA14+SA_salmon/PA14+PAO1+BC+SA/aux_info/

seqkit grep -n -f modified_names.txt /scratch/gpfs/ziyangc/PASA/ScRNAseq/new/rawdata/condition_demultiplex_full/PA14+SA_salmon/PA14+PAO1+BC/aux_info/unmapped_reads.fastq.gz \
-o unmapped_reads.fastq

gzip unmapped_reads.fastq