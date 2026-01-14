# M3-seq Bacterial scRNA-seq Pipeline

A SLURM-based pipeline for processing bacterial single-cell RNA sequencing data with combinatorial indexing (M3-seq).

## Overview

This pipeline processes bacterial scRNA-seq data through the following steps:

1. **Preprocessing** - Trim reads to extract UMI and cell barcodes
2. **Demultiplexing** - Sequential demultiplexing by i7 index, then plate barcode
3. **Quantification** - Alignment and UMI counting with kallisto-bustools or alevin-fry

## Read Structure

The pipeline expects 4 FASTQ files from Illumina sequencing:

| Read | Length | Structure | Content |
|------|--------|-----------|---------|
| Read 1 | 26 bp | `[8 UMI][13 plate BC][5 skip]` | UMI (rx tag) + plate barcode for demux |
| Read 2 | 8 bp | `[8 i7]` | i7 index barcode |
| Read 3 | 30 bp | `[14 skip][16 cell BC]` | 10x cell barcode (r2 tag, reverse-complemented) |
| Read 4 | 74 bp | `[74 transcript]` | Transcript sequence |

Concatenated structure: `8M13B5S8B14S16M74T`

## Pipeline Steps

```
┌─────────────────────────────────────────────────────────────────┐
│ 01_preprocess.slurm                                             │
│   • Trim R1 → extract UMI (8bp) and plate barcode region        │
│   • Trim R3 → extract cell barcode (16bp)                       │
│   • Combine cell barcode + UMI → 24bp barcode read              │
└─────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────┐
│ 02_demultiplex.slurm                                            │
│   • Generate barcode references from sample sheet               │
│   • Demux by i7 (Read 2) → split all reads by i7 index          │
│   • Demux by plate barcode (Read 1) → split by well             │
│   • Merge samples by LIBRARY_NAME for quantification            │
└─────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────┐
│ 03_quantify.slurm                                               │
│   • Optional: fastp quality control                             │
│   • Quantify with kallisto-bustools OR alevin-fry               │
│   • Generate per-library statistics                             │
└─────────────────────────────────────────────────────────────────┘
```

## Requirements

### Software Dependencies

- **cutadapt** (≥4.0) - Read trimming and demultiplexing
- **kallisto** + **bustools** (kb-python) - Quantification option 1
- **simpleaf** (alevin-fry) - Quantification option 2
- **fastp** (optional) - Quality control
- **pigz** (optional) - Parallel gzip compression
- **Python 3** with packages: `pyyaml`, `biopython`

### Conda Environment Setup

```bash
conda create -n scseq python=3.10
conda activate scseq
conda install -c bioconda cutadapt kb-python fastp pigz
pip install pyyaml biopython

# For alevin-fry (optional)
conda install -c bioconda simpleaf
```

## Directory Structure

```
pipeline/
├── config.yaml              # Configuration file (edit this)
├── samples.csv              # Sample sheet
├── run_pipeline.sh          # Master submission script
├── scripts/
│   ├── 01_preprocess.slurm
│   ├── 02_demultiplex.slurm
│   └── 03_quantify.slurm
└── utils/
    └── combine_fastq.py
```

## Configuration

Edit `config.yaml` before running. Key sections:

### SLURM Settings

```yaml
slurm:
  partition: "your_partition"
  account: "your_account"       # Leave empty if not required

environment:
  module_load:
    - "anaconda3/2024.10"       # Adjust for your cluster
  conda_env: "scseq"
```

### Input/Output Paths

```yaml
paths:
  raw_data_dir: "/path/to/rawdata"
  read1: "Sample_Read_1.fastq.gz"
  read2: "Sample_Read_2.fastq.gz"    # i7 index
  read3: "Sample_Read_3.fastq.gz"
  read4: "Sample_Read_4.fastq.gz"    # transcript

  sample_sheet: "samples.csv"
  output_dir: "/path/to/output"

  whitelist: "/path/to/rc_737K-cratac-v1.txt"  # 10x whitelist (RC)
  genome_dir: "/path/to/genome"
  kallisto_index: "kallisto/index.idx"
  kallisto_t2g: "kallisto/t2g.txt"
```

### Parallelization

```yaml
resources:
  preprocess:
    cpus: 4
    threads: 4

  demultiplex:
    cpus: 16
    i7_threads: 8
    plate_threads: 4
    pigz_threads: 4

  quantify:
    cpus: 8
    threads: 8
```

### Demultiplexing Settings

```yaml
demultiplex:
  i7_error_rate: 0           # Exact match for i7
  plate_error_rate: 0.1      # ~1 mismatch allowed for 13bp plate BC
```

### Quantification Tool

```yaml
quantification:
  tool: "kallisto"           # Options: "kallisto" or "alevin"
```

## Sample Sheet Format

The sample sheet (`samples.csv`) maps barcode combinations to samples:

| Column | Description |
|--------|-------------|
| SKIP | 0 = process, 1 = skip |
| LIBRARY_NAME | Experimental condition (samples merged by this) |
| SAMPLE_NAME | Unique sample identifier (e.g., well position) |
| BARCODE_1 | 13bp plate barcode (from Read 1) |
| BARCODE_2 | 8bp i7 index barcode (from Read 2) |

Example:
```csv
SKIP,LIBRARY_NAME,SAMPLE_NAME,BARCODE_1,BARCODE_2
0,Lib1_MG1655_Exponential,Lib1_MG1655_Exponential_A1,AAGTGATTAGCAA,CGTACTAG
0,Lib1_MG1655_Exponential,Lib1_MG1655_Exponential_A2,AGAATCCCCCTAA,CGTACTAG
0,Lib1_MG1655_Tet,Lib1_MG1655_Tet_A9,AACGTTCTGTCGA,CGTACTAG
```

## Usage

### 1. Prepare Reference Files

Build kallisto or alevin index:

```bash
# Kallisto
kb ref -i index.idx -g t2g.txt -f1 transcriptome.fa genome.fna genes.gtf

# Alevin-fry
simpleaf index --output genome_dir --ref-seq transcriptome.fa
```

### 2. Configure Pipeline

```bash
cd pipeline
cp config.yaml my_config.yaml
# Edit my_config.yaml with your paths
```

### 3. Run Pipeline

```bash
./run_pipeline.sh my_config.yaml
```

### 4. Monitor Jobs

```bash
squeue -u $USER
```

### 5. Check Logs

```bash
ls /path/to/output/logs/
tail -f /path/to/output/logs/01_preprocess_*.out
```

## Output Structure

```
output/
├── 01_preprocessed/
│   ├── trimmed_read1_umi.fastq.gz
│   ├── trimmed_read1_platebc.fastq.gz
│   ├── trimmed_read3_cellbc.fastq.gz
│   ├── combined_barcode_umi.fastq.gz    # Cell BC + UMI (24bp)
│   ├── read2_i7.fastq.gz
│   └── read4_transcript.fastq.gz
│
├── 02_demultiplexed/
│   ├── barcodes/                         # Generated barcode files
│   ├── i7_demux/                         # Reads split by i7
│   └── plate_demux/                      # Reads split by plate BC
│
├── 03_merged_libraries/
│   ├── {LIBRARY_NAME}.barcode_umi.fastq.gz
│   └── {LIBRARY_NAME}.transcript.fastq.gz
│
├── 04_qc/                                # If QC enabled
│   ├── {LIBRARY_NAME}.fastp.html
│   └── {LIBRARY_NAME}.fastp.json
│
├── 05_quantification/
│   └── {LIBRARY_NAME}/
│       ├── counts_unfiltered/
│       ├── counts_filtered/
│       └── *.h5ad                        # AnnData output
│
├── 06_stats/
│   ├── pipeline_stats.tsv
│   └── pipeline_report.txt
│
└── logs/
    ├── 01_preprocess_*.out
    ├── 02_demultiplex_*.out
    └── 03_quantify_*.out
```

## Troubleshooting

### Job fails immediately

- Check SLURM partition name and account in config
- Verify module names match your cluster

### No reads after demultiplexing

- Verify barcode sequences in sample sheet match your library prep
- Check `cutadapt_*.log` files for barcode matching statistics
- Try increasing `plate_error_rate` if barcodes have synthesis errors

### Low mapping rates

- Verify genome/transcriptome reference matches your organism
- Check strandedness setting in quantification config
- Review fastp QC reports if enabled

### Memory errors

- Increase `mem` in resources section
- Reduce `threads` to lower memory footprint

### Permission errors

```bash
chmod +x run_pipeline.sh
chmod +x utils/combine_fastq.py
```

## Citation

If you use this pipeline, please cite the relevant tools:

- cutadapt: Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10-12.
- kallisto/bustools: Melsted, P., et al. (2021). Modular, efficient and constant-memory single-cell RNA-seq preprocessing. Nature Biotechnology.
- alevin-fry: He, D., et al. (2022). Alevin-fry unlocks rapid, accurate and memory-frugal quantification of single-cell RNA-seq data. Nature Methods.

## License

See LICENSE file for details.
