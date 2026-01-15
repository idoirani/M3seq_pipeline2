#!/bin/bash
# =============================================================================
# 02_demultiplex.slurm - Demultiplex by i7 and plate barcodes
# =============================================================================
# This script:
#   1. Generates barcode reference files from sample sheet
#   2. Demultiplexes by i7 barcode (Read 2) using cutadapt
#   3. Demultiplexes by plate barcode (Read 1) using cutadapt
#   4. Merges samples by LIBRARY_NAME for downstream quantification
# =============================================================================

#SBATCH --job-name=m3seq_demultiplex
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=12:00:00

# Note: Actual resource values are overridden by run_pipeline.sh from config.yaml

set -euo pipefail

# -----------------------------------------------------------------------------
# Load configuration
# -----------------------------------------------------------------------------
echo "=============================================="
echo "M3-seq Demultiplexing"
echo "=============================================="
echo "Start time: $(date)"
echo "Config file: $CONFIG_FILE"
echo ""

# Simple YAML parser function
parse_yaml() {
    python3 << EOF
import yaml
with open("$CONFIG_FILE", 'r') as f:
    config = yaml.safe_load(f)
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

MODULES=$(parse_yaml "environment.module_load")
for mod in $MODULES; do
    module load "$mod"
done

CONDA_ENV=$(parse_yaml "environment.conda_env")
set +u  # Disable unbound variable check for conda activation
conda activate "$CONDA_ENV"
set -u  # Re-enable

# -----------------------------------------------------------------------------
# Get paths and settings from config
# -----------------------------------------------------------------------------
OUTPUT_DIR=$(parse_yaml "paths.output_dir")
SAMPLE_SHEET=$(parse_yaml "paths.sample_sheet")

# Handle relative sample sheet path
if [[ ! "$SAMPLE_SHEET" = /* ]]; then
    SAMPLE_SHEET="${PIPELINE_DIR}/${SAMPLE_SHEET}"
fi

I7_ERROR_RATE=$(parse_yaml "demultiplex.i7_error_rate")
PLATE_ERROR_RATE=$(parse_yaml "demultiplex.plate_error_rate")

# Parallelization settings
I7_THREADS=$(parse_yaml "resources.demultiplex.i7_threads")
PLATE_THREADS=$(parse_yaml "resources.demultiplex.plate_threads")
PIGZ_THREADS=$(parse_yaml "resources.demultiplex.pigz_threads")

PREPROCESS_DIR="${OUTPUT_DIR}/01_preprocessed"
DEMUX_DIR="${OUTPUT_DIR}/02_demultiplexed"
BARCODES_DIR="${DEMUX_DIR}/barcodes"
I7_DEMUX_DIR="${DEMUX_DIR}/i7_demux"
PLATE_DEMUX_DIR="${DEMUX_DIR}/plate_demux"
MERGED_DIR="${OUTPUT_DIR}/03_merged_libraries"

mkdir -p "$BARCODES_DIR" "$I7_DEMUX_DIR" "$PLATE_DEMUX_DIR" "$MERGED_DIR"

echo "Sample sheet: $SAMPLE_SHEET"
echo "Output directory: $DEMUX_DIR"
echo ""

# -----------------------------------------------------------------------------
# Step 1: Generate barcode reference files from sample sheet
# -----------------------------------------------------------------------------
echo "=============================================="
echo "Step 1: Generating barcode reference files"
echo "=============================================="

python3 << 'PYTHON_SCRIPT'
import csv
import os

sample_sheet = os.environ.get('SAMPLE_SHEET', '')
barcodes_dir = os.environ.get('BARCODES_DIR', '')

# Read sample sheet
samples = []
with open(sample_sheet, 'r', encoding='utf-8-sig') as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row['SKIP'] == '0':
            samples.append(row)

print(f"Loaded {len(samples)} samples (excluding skipped)")

# Get unique i7 barcodes
i7_barcodes = {}
for s in samples:
    bc = s['BARCODE_2']
    if bc not in i7_barcodes:
        i7_barcodes[bc] = f"i7_{bc}"

# Write i7 barcode FASTA
i7_fasta = os.path.join(barcodes_dir, 'i7_barcodes.fasta')
with open(i7_fasta, 'w') as f:
    for bc, name in i7_barcodes.items():
        f.write(f">{name}\n^{bc}\n")
print(f"Created i7 barcode file: {i7_fasta} ({len(i7_barcodes)} barcodes)")

# Get unique plate barcodes (BARCODE_1) per i7
# Group samples by i7
samples_by_i7 = {}
for s in samples:
    i7 = s['BARCODE_2']
    if i7 not in samples_by_i7:
        samples_by_i7[i7] = []
    samples_by_i7[i7].append(s)

# Write plate barcode FASTA for each i7 (or combined)
plate_fasta = os.path.join(barcodes_dir, 'plate_barcodes.fasta')
with open(plate_fasta, 'w') as f:
    seen = set()
    for s in samples:
        bc = s['BARCODE_1']
        name = s['SAMPLE_NAME']
        if (bc, name) not in seen:
            # Use 3' anchored adapter for plate barcode matching
            f.write(f">{name}\n{bc}$\n")
            seen.add((bc, name))
print(f"Created plate barcode file: {plate_fasta} ({len(seen)} entries)")

# Write sample-to-library mapping for merging step
mapping_file = os.path.join(barcodes_dir, 'sample_library_mapping.tsv')
with open(mapping_file, 'w') as f:
    f.write("SAMPLE_NAME\tLIBRARY_NAME\tBARCODE_1\tBARCODE_2\n")
    for s in samples:
        f.write(f"{s['SAMPLE_NAME']}\t{s['LIBRARY_NAME']}\t{s['BARCODE_1']}\t{s['BARCODE_2']}\n")
print(f"Created mapping file: {mapping_file}")

# Write list of unique libraries
libraries = sorted(set(s['LIBRARY_NAME'] for s in samples))
lib_file = os.path.join(barcodes_dir, 'libraries.txt')
with open(lib_file, 'w') as f:
    for lib in libraries:
        f.write(f"{lib}\n")
print(f"Created library list: {lib_file} ({len(libraries)} libraries)")

# Write i7 list for iteration
i7_file = os.path.join(barcodes_dir, 'i7_list.txt')
with open(i7_file, 'w') as f:
    for bc in i7_barcodes.keys():
        f.write(f"{bc}\n")
print(f"Created i7 list: {i7_file}")
PYTHON_SCRIPT

export SAMPLE_SHEET BARCODES_DIR

echo ""

# -----------------------------------------------------------------------------
# Step 2: Demultiplex by i7 barcode
# -----------------------------------------------------------------------------
echo "=============================================="
echo "Step 2: Demultiplexing by i7 barcode"
echo "=============================================="

# Input files from preprocessing
BARCODE_UMI="${PREPROCESS_DIR}/combined_barcode_umi.fastq.gz"
TRANSCRIPT="${PREPROCESS_DIR}/read4_transcript.fastq.gz"
I7_READ="${PREPROCESS_DIR}/read2_i7.fastq.gz"
PLATE_BC_READ="${PREPROCESS_DIR}/trimmed_read1_platebc.fastq.gz"

echo "Using i7 read for demultiplexing: $I7_READ"
echo "Paired files to demultiplex:"
echo "  - Barcode+UMI: $BARCODE_UMI"
echo "  - Transcript: $TRANSCRIPT"
echo "  - Plate BC read: $PLATE_BC_READ"
echo ""

# Demultiplex by i7 - this splits all reads by i7 barcode
# We need to demux 3 files in sync: barcode_umi, transcript, and plate_bc
# cutadapt handles 2 files at a time (paired), so we run parallel passes

echo "Running cutadapt for i7 demultiplexing..."
echo "Using $I7_THREADS threads per i7 demux job"
echo ""

mkdir -p "${I7_DEMUX_DIR}/barcode_umi" "${I7_DEMUX_DIR}/transcript" "${I7_DEMUX_DIR}/plate_bc"
mkdir -p "${I7_DEMUX_DIR}/i7_pass1" "${I7_DEMUX_DIR}/i7_pass2" "${I7_DEMUX_DIR}/i7_pass3"

# Pass 1: demux i7 + barcode_umi pair
cutadapt -e "${I7_ERROR_RATE}" -j "${I7_THREADS}" --no-indels \
    -g "file:${BARCODES_DIR}/i7_barcodes.fasta" \
    -o "${I7_DEMUX_DIR}/i7_pass1/{name}.fastq.gz" \
    -p "${I7_DEMUX_DIR}/barcode_umi/{name}.fastq.gz" \
    "$I7_READ" "$BARCODE_UMI" \
    2>&1 | tee "${DEMUX_DIR}/cutadapt_i7_pass1.log" &
PID1=$!

# Pass 2: demux i7 + transcript pair
cutadapt -e "${I7_ERROR_RATE}" -j "${I7_THREADS}" --no-indels \
    -g "file:${BARCODES_DIR}/i7_barcodes.fasta" \
    -o "${I7_DEMUX_DIR}/i7_pass2/{name}.fastq.gz" \
    -p "${I7_DEMUX_DIR}/transcript/{name}.fastq.gz" \
    "$I7_READ" "$TRANSCRIPT" \
    2>&1 | tee "${DEMUX_DIR}/cutadapt_i7_pass2.log" &
PID2=$!

# Pass 3: demux i7 + plate_bc pair
cutadapt -e "${I7_ERROR_RATE}" -j "${I7_THREADS}" --no-indels \
    -g "file:${BARCODES_DIR}/i7_barcodes.fasta" \
    -o "${I7_DEMUX_DIR}/i7_pass3/{name}.fastq.gz" \
    -p "${I7_DEMUX_DIR}/plate_bc/{name}.fastq.gz" \
    "$I7_READ" "$PLATE_BC_READ" \
    2>&1 | tee "${DEMUX_DIR}/cutadapt_i7_pass3.log" &
PID3=$!

wait $PID1 $PID2 $PID3

# Clean up intermediate i7 files (we only need the paired outputs)
rm -rf "${I7_DEMUX_DIR}/i7_pass1" "${I7_DEMUX_DIR}/i7_pass2" "${I7_DEMUX_DIR}/i7_pass3"

echo "i7 demultiplexing complete."
echo ""

# -----------------------------------------------------------------------------
# Step 3: Demultiplex by plate barcode (within each i7 group)
# -----------------------------------------------------------------------------
echo "=============================================="
echo "Step 3: Demultiplexing by plate barcode"
echo "=============================================="

# Read i7 list
while IFS= read -r I7_BC; do
    I7_NAME="i7_${I7_BC}"
    echo "Processing i7 group: $I7_NAME"

    # Check if files exist for this i7
    if [[ ! -f "${I7_DEMUX_DIR}/barcode_umi/${I7_NAME}.fastq.gz" ]]; then
        echo "  Skipping - no reads for this i7"
        continue
    fi

    # Create output directories for this i7
    mkdir -p "${PLATE_DEMUX_DIR}/${I7_NAME}/barcode_umi"
    mkdir -p "${PLATE_DEMUX_DIR}/${I7_NAME}/transcript"

    # Demultiplex by plate barcode - run two passes to keep barcode_umi and transcript in sync
    mkdir -p "${PLATE_DEMUX_DIR}/${I7_NAME}/plate_pass1" "${PLATE_DEMUX_DIR}/${I7_NAME}/plate_pass2"

    # Pass 1: plate_bc + barcode_umi
    cutadapt -e "${PLATE_ERROR_RATE}" -j "${PLATE_THREADS}" \
        -a "file:${BARCODES_DIR}/plate_barcodes.fasta" \
        -o "${PLATE_DEMUX_DIR}/${I7_NAME}/plate_pass1/{name}.fastq.gz" \
        -p "${PLATE_DEMUX_DIR}/${I7_NAME}/barcode_umi/{name}.fastq.gz" \
        "${I7_DEMUX_DIR}/plate_bc/${I7_NAME}.fastq.gz" \
        "${I7_DEMUX_DIR}/barcode_umi/${I7_NAME}.fastq.gz" \
        2>&1 | tee "${DEMUX_DIR}/cutadapt_plate_${I7_NAME}_pass1.log" &
    PID_PLATE1=$!

    # Pass 2: plate_bc + transcript
    cutadapt -e "${PLATE_ERROR_RATE}" -j "${PLATE_THREADS}" \
        -a "file:${BARCODES_DIR}/plate_barcodes.fasta" \
        -o "${PLATE_DEMUX_DIR}/${I7_NAME}/plate_pass2/{name}.fastq.gz" \
        -p "${PLATE_DEMUX_DIR}/${I7_NAME}/transcript/{name}.fastq.gz" \
        "${I7_DEMUX_DIR}/plate_bc/${I7_NAME}.fastq.gz" \
        "${I7_DEMUX_DIR}/transcript/${I7_NAME}.fastq.gz" \
        2>&1 | tee "${DEMUX_DIR}/cutadapt_plate_${I7_NAME}_pass2.log" &
    PID_PLATE2=$!

    # Wait for both passes to complete
    wait $PID_PLATE1 $PID_PLATE2

    # Clean up intermediate files
    rm -rf "${PLATE_DEMUX_DIR}/${I7_NAME}/plate_pass1" "${PLATE_DEMUX_DIR}/${I7_NAME}/plate_pass2"

    echo "  Plate demultiplexing complete for $I7_NAME"
done < "${BARCODES_DIR}/i7_list.txt"

echo ""

# -----------------------------------------------------------------------------
# Step 4: Merge samples by LIBRARY_NAME
# -----------------------------------------------------------------------------
echo "=============================================="
echo "Step 4: Merging samples by LIBRARY_NAME"
echo "=============================================="

# Read library list and merge
while IFS= read -r LIBRARY; do
    echo "Merging library: $LIBRARY"

    # Find all samples belonging to this library
    SAMPLE_LIST=$(grep -P "\t${LIBRARY}\t" "${BARCODES_DIR}/sample_library_mapping.tsv" | cut -f1)

    # Initialize merged files
    MERGED_BC="${MERGED_DIR}/${LIBRARY}.barcode_umi.fastq.gz"
    MERGED_TX="${MERGED_DIR}/${LIBRARY}.transcript.fastq.gz"

    # Collect files to merge
    BC_FILES=""
    TX_FILES=""

    for SAMPLE in $SAMPLE_LIST; do
        # Get the i7 barcode for this sample
        I7_BC=$(grep -P "^${SAMPLE}\t" "${BARCODES_DIR}/sample_library_mapping.tsv" | cut -f4)
        I7_NAME="i7_${I7_BC}"

        BC_FILE="${PLATE_DEMUX_DIR}/${I7_NAME}/barcode_umi/${SAMPLE}.fastq.gz"
        TX_FILE="${PLATE_DEMUX_DIR}/${I7_NAME}/transcript/${SAMPLE}.fastq.gz"

        if [[ -f "$BC_FILE" ]]; then
            BC_FILES="$BC_FILES $BC_FILE"
            TX_FILES="$TX_FILES $TX_FILE"
        else
            echo "  Warning: File not found for sample $SAMPLE"
        fi
    done

    # Merge files
    if [[ -n "$BC_FILES" ]]; then
        echo "  Concatenating $(echo $BC_FILES | wc -w) sample files..."
        # Use pigz for parallel compression if available, fallback to gzip
        if command -v pigz &> /dev/null; then
            zcat $BC_FILES | pigz -p "${PIGZ_THREADS}" > "$MERGED_BC" &
            PID_BC=$!
            zcat $TX_FILES | pigz -p "${PIGZ_THREADS}" > "$MERGED_TX" &
            PID_TX=$!
            wait $PID_BC $PID_TX
        else
            zcat $BC_FILES | gzip > "$MERGED_BC"
            zcat $TX_FILES | gzip > "$MERGED_TX"
        fi

        # Count reads
        NREADS=$(zcat "$MERGED_BC" | wc -l)
        NREADS=$((NREADS / 4))
        echo "  Merged reads: $NREADS"
    else
        echo "  Warning: No files found for library $LIBRARY"
    fi
done < "${BARCODES_DIR}/libraries.txt"

echo ""

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo "=============================================="
echo "Demultiplexing Summary"
echo "=============================================="
echo ""
echo "Output directories:"
echo "  - i7 demultiplexed: $I7_DEMUX_DIR"
echo "  - Plate demultiplexed: $PLATE_DEMUX_DIR"
echo "  - Merged by library: $MERGED_DIR"
echo ""

echo "Merged libraries:"
for f in "${MERGED_DIR}"/*.barcode_umi.fastq.gz; do
    if [[ -f "$f" ]]; then
        LIB=$(basename "$f" .barcode_umi.fastq.gz)
        NREADS=$(zcat "$f" | wc -l)
        NREADS=$((NREADS / 4))
        echo "  $LIB: $NREADS reads"
    fi
done

echo ""
echo "End time: $(date)"
echo "Demultiplexing complete!"
