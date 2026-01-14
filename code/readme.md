read1 has 5'-UMI[8]-PlateBarcode[13]-Trash[5]-3'
read2 has i5
read3 has 5'-Trash[14]-CellBarcode(reverse compliment)[16]-3'
read4 has sequence


1. wget data
2. md5_checksum.sh -> md5check.ipynb
3. trim.slurm
   1. remove the 3' 18 nt from read1, so trimmed_read1.fastq.gz is UMI
   2. remove the 5' 14 nt from read3, so trimmed_read3.fastq.gz is cellbarcode from 10x
   3. remove th3 5' 8 nt from read1, so trimmed_read1_platebarcode.fastq.gz is platebarcode+trash[5]
4. combine_fast.slurm
   1. combine trimmed_read1 and trimmed_read3 to get the format 5'-cellbarcode[16]-UMI[8]-3' named trimmed_read3_and_1.fastq.gz
5. whitelist.slurm (can skip, use rc-737K instead)
   1. generate cellbarcode whitelist using trimmed_read3_and_1.fastq.gz
6. barcode_for_demult.ipynb
   1. fasta format file for well demultiplexing
7. demultiplex.slurm
   1. demultiplex read3_and_1 and read4 using well_demultiple file
8. concatenate.ipynb -> concatenate.slurm
   1. concatenate well_demultiplexed files using condition
   2. the output from demultiplex.slurm can be removed after this is done
9.  