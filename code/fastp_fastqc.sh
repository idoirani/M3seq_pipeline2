fastp -i in.fastq.gz -I in.fastq.gz -o out.fastq.gz -O out.fastq.gz --trim_poly_g --detect_adapter_for_pe
fastqc ./PA14+SA.4.fastp.fastq.gz -o ../fastqc/PA14+SA_fastp

fastp -i in.fastq.gz -o out.fastq.gz --trim_ploy_g --json ./fastp/report.json --html ./fastp/report.html