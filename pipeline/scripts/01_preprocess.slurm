#!/bin/bash
# =============================================================================
# 01_preprocess.slurm - Trim reads and combine barcodes
# =============================================================================
# This script:
#   1. Trims Read 1 to extract UMI (first 8bp)
#   2. Trims Read 1 to extract plate barcode region (for demultiplexing)
#   3. Trims Read 3 to extract cell barcode (skip first 14bp, keep 16bp)
#   4. Combines cell barcode + UMI into single read for downstream processing
# =============================================================================

#SBATCH --job-name=m3seq_preprocess
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00

# Note: Actual resource values are overridden by run_pipeline.sh from config.yaml

set -eo pipefail

# Set PS1 for non-interactive shell (required for conda)
export PS1=${PS1:-""}

# -----------------------------------------------------------------------------
# Load configuration
# -----------------------------------------------------------------------------
echo "=============================================="
echo "M3-seq Preprocessing"
echo "=============================================="
echo "Start time: $(date)"

# Validate required environment variables
if [[ -z "${CONFIG_FILE:-}" ]]; then
    echo "ERROR: CONFIG_FILE environment variable not set"
    exit 1
fi
if [[ -z "${PIPELINE_DIR:-}" ]]; then
    echo "ERROR: PIPELINE_DIR environment variable not set"
    exit 1
fi
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

echo "Config file: $CONFIG_FILE"
echo "Pipeline dir: $PIPELINE_DIR"
echo ""

# Simple YAML parser function
parse_yaml() {
    python3 << EOF
import yaml
import sys

with open("$CONFIG_FILE", 'r') as f:
    config = yaml.safe_load(f)

# Navigate nested keys
keys = "$1".split('.')
value = config
for key in keys:
    value = value[key]

if isinstance(value, list):
    print(' '.join(str(v) for v in value))
else:
    print(value)
EOF
}

# -----------------------------------------------------------------------------
# Load environment
# -----------------------------------------------------------------------------
module purge

# Load modules from config
MODULES=$(parse_yaml "environment.module_load")
echo "Modules to load: $MODULES"
for mod in $MODULES; do
    echo "Loading module: $mod"
    module load "$mod"
done

CONDA_ENV=$(parse_yaml "environment.conda_env")
echo "Activating conda environment: $CONDA_ENV"
source activate "$CONDA_ENV" || conda activate "$CONDA_ENV"

# Verify key commands are available
echo "Checking required commands..."
command -v cutadapt >/dev/null 2>&1 || { echo "ERROR: cutadapt not found"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found"; exit 1; }
echo "All required commands found."

# -----------------------------------------------------------------------------
# Get paths from config
# -----------------------------------------------------------------------------
RAW_DATA_DIR=$(parse_yaml "paths.raw_data_dir")
READ1_FILE=$(parse_yaml "paths.read1")
READ2_FILE=$(parse_yaml "paths.read2")
READ3_FILE=$(parse_yaml "paths.read3")
READ4_FILE=$(parse_yaml "paths.read4")
OUTPUT_DIR=$(parse_yaml "paths.output_dir")

# Read structure
UMI_LEN=$(parse_yaml "read_structure.read1_umi_len")
PLATE_BC_START=$(parse_yaml "read_structure.read1_plate_bc_start")
PLATE_BC_LEN=$(parse_yaml "read_structure.read1_plate_bc_len")
READ1_SKIP_LEN=$(parse_yaml "read_structure.read1_skip_len")
READ3_SKIP_LEN=$(parse_yaml "read_structure.read3_skip_len")
CELL_BC_LEN=$(parse_yaml "read_structure.read3_cell_bc_len")

# Parallelization settings
THREADS=$(parse_yaml "resources.preprocess.threads")

# Full paths
READ1="${RAW_DATA_DIR}/${READ1_FILE}"
READ2="${RAW_DATA_DIR}/${READ2_FILE}"
READ3="${RAW_DATA_DIR}/${READ3_FILE}"
READ4="${RAW_DATA_DIR}/${READ4_FILE}"

# Create output directories
PREPROCESS_DIR="${OUTPUT_DIR}/01_preprocessed"
mkdir -p "$PREPROCESS_DIR"

echo ""
echo "Input files:"
echo "  Read 1: $READ1"
echo "  Read 2: $READ2"
echo "  Read 3: $READ3"
echo "  Read 4: $READ4"
echo ""
echo "Output directory: $PREPROCESS_DIR"
echo ""

# Validate input files
for f in "$READ1" "$READ2" "$READ3" "$READ4"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Input file not found: $f"
        exit 1
    fi
done

# -----------------------------------------------------------------------------
# Step 1: Trim reads
# -----------------------------------------------------------------------------
echo "=============================================="
echo "Step 1: Trimming reads"
echo "=============================================="

# Calculate trim lengths for Read 1
# Read 1 structure: [8 UMI][13 plate BC][5 skip] = 26bp total
# For UMI extraction: keep first 8bp, trim last 18bp
TRIM_R1_END=$((PLATE_BC_LEN + READ1_SKIP_LEN))  # 13 + 5 = 18

# For plate barcode extraction (for demux): skip first 8bp (UMI)
# This gives us: [13 plate BC][5 skip] which cutadapt will use for demux

echo "Using $THREADS threads for trimming"
echo ""

echo "Trimming Read 1 for UMI extraction (keep first ${UMI_LEN}bp)..."
cutadapt -j ${THREADS} -u -${TRIM_R1_END} \
    -o "${PREPROCESS_DIR}/trimmed_read1_umi.fastq.gz" \
    "$READ1" &
PID1=$!

echo "Trimming Read 1 for plate barcode (skip first ${UMI_LEN}bp)..."
cutadapt -j ${THREADS} -u ${UMI_LEN} \
    -o "${PREPROCESS_DIR}/trimmed_read1_platebc.fastq.gz" \
    "$READ1" &
PID2=$!

echo "Trimming Read 3 for cell barcode (skip first ${READ3_SKIP_LEN}bp)..."
cutadapt -j ${THREADS} -u ${READ3_SKIP_LEN} \
    -o "${PREPROCESS_DIR}/trimmed_read3_cellbc.fastq.gz" \
    "$READ3" &
PID3=$!

# Copy Read 2 (i7) and Read 4 (transcript) - no trimming needed
echo "Linking Read 2 (i7 index)..."
ln -sf "$READ2" "${PREPROCESS_DIR}/read2_i7.fastq.gz" &
PID4=$!

echo "Linking Read 4 (transcript)..."
ln -sf "$READ4" "${PREPROCESS_DIR}/read4_transcript.fastq.gz" &
PID5=$!

# Wait for all trimming jobs
wait $PID1 $PID2 $PID3 $PID4 $PID5
echo "Trimming complete."
echo ""

# -----------------------------------------------------------------------------
# Step 2: Combine cell barcode + UMI
# -----------------------------------------------------------------------------
echo "=============================================="
echo "Step 2: Combining cell barcode + UMI"
echo "=============================================="

# Combine: [16bp cell barcode][8bp UMI] = 24bp barcode read
# This format is compatible with kallisto-bustools and alevin-fry

echo "Combining trimmed_read3_cellbc + trimmed_read1_umi..."
python3 "${PIPELINE_DIR}/utils/combine_fastq.py" \
    "${PREPROCESS_DIR}/trimmed_read3_cellbc.fastq.gz" \
    "${PREPROCESS_DIR}/trimmed_read1_umi.fastq.gz" \
    "${PREPROCESS_DIR}/combined_barcode_umi.fastq.gz"

echo "Combined barcode file created: ${PREPROCESS_DIR}/combined_barcode_umi.fastq.gz"
echo ""

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo "=============================================="
echo "Preprocessing Summary"
echo "=============================================="
echo ""
echo "Output files created:"
echo "  - trimmed_read1_umi.fastq.gz      : UMI only (${UMI_LEN}bp)"
echo "  - trimmed_read1_platebc.fastq.gz  : Plate barcode region for demux"
echo "  - trimmed_read3_cellbc.fastq.gz   : Cell barcode (${CELL_BC_LEN}bp)"
echo "  - read2_i7.fastq.gz               : i7 index (link)"
echo "  - read4_transcript.fastq.gz       : Transcript (link)"
echo "  - combined_barcode_umi.fastq.gz   : Cell BC + UMI ($((CELL_BC_LEN + UMI_LEN))bp)"
echo ""

# Count reads
echo "Read counts:"
NREADS=$(zcat "${PREPROCESS_DIR}/combined_barcode_umi.fastq.gz" | wc -l)
NREADS=$((NREADS / 4))
echo "  Total reads: $NREADS"
echo ""

echo "End time: $(date)"
echo "Preprocessing complete!"
