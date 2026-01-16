#!/bin/bash
# =============================================================================
# Dump first 100 lines of each FASTQ file for debugging
# =============================================================================
# Usage: ./dump_read_samples.sh
# Output files will be created in the current directory
# =============================================================================

# Raw data directory
RAW_DIR="/scratch/gpfs/GITAI/ido/M3seq/raw_data/bw_scifi09"
PREPROCESS_DIR="/scratch/gpfs/GITAI/ido/M3seq/output_scifi09_pipe2/01_preprocessed"

# Output directory (current directory or specify one)
OUT_DIR="${1:-.}"

echo "Dumping first 100 lines from each FASTQ file..."
echo "Output directory: $OUT_DIR"
echo ""

# Read 1: 26bp (8bp UMI + 13bp plate barcode + 5bp skip)
echo "Extracting Read 1 (raw)..."
zcat "${RAW_DIR}/read1.fastq.gz" | head -100 > "${OUT_DIR}/read1_sample.txt"

# Read 2: 8bp i7 index
echo "Extracting Read 2 (raw)..."
zcat "${RAW_DIR}/read2.fastq.gz" | head -100 > "${OUT_DIR}/read2_sample.txt"

# Read 3: 30bp (14bp skip + 16bp cell barcode)
echo "Extracting Read 3 (raw)..."
zcat "${RAW_DIR}/read3.fastq.gz" | head -100 > "${OUT_DIR}/read3_sample.txt"

# Read 4: 74bp transcript
echo "Extracting Read 4 (raw)..."
zcat "${RAW_DIR}/read4.fastq.gz" | head -100 > "${OUT_DIR}/read4_sample.txt"

# Trimmed Read 1 for plate barcode (after UMI removal)
echo "Extracting trimmed Read 1 plate BC..."
zcat "${PREPROCESS_DIR}/trimmed_read1_platebc.fastq.gz" | head -100 > "${OUT_DIR}/trimmed_read1_platebc_sample.txt"

echo ""
echo "Done! Created files:"
ls -la "${OUT_DIR}"/*.txt

echo ""
echo "To copy to local machine:"
echo "  scp della:${OUT_DIR}/*_sample.txt ."
