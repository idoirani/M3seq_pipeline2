#!/bin/bash
# =============================================================================
# M3-seq Bacterial scRNA-seq Pipeline - Master Orchestration Script
# =============================================================================
# Usage: ./run_pipeline.sh [config.yaml]
#
# This script submits all pipeline jobs with proper dependencies.
# Jobs will run sequentially: preprocess -> demultiplex -> quantify
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Parse arguments and setup
# -----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${1:-${SCRIPT_DIR}/config.yaml}"

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    echo "Usage: $0 [config.yaml]"
    exit 1
fi

echo "=============================================="
echo "M3-seq Pipeline"
echo "=============================================="
echo "Config: $CONFIG_FILE"
echo "Script directory: $SCRIPT_DIR"
echo ""

# -----------------------------------------------------------------------------
# Parse YAML config (simple bash parsing for key values)
# -----------------------------------------------------------------------------
parse_yaml() {
    local yaml_file="$1"
    local key="$2"
    grep -E "^\s*${key}:" "$yaml_file" | head -1 | sed 's/.*:\s*//' | sed 's/ *#.*$//' | tr -d '"' | tr -d "'"
}

parse_yaml_nested() {
    local yaml_file="$1"
    local section="$2"
    local key="$3"
    awk -v section="$section" -v key="$key" '
        $0 ~ "^"section":" { in_section=1; next }
        in_section && /^[a-z]/ { in_section=0 }
        in_section && $0 ~ "^  "key":" {
            gsub(/.*: */, "");
            gsub(/["\047]/, "");
            gsub(/ *#.*$/, "");  # Remove comments
            print;
            exit
        }
    ' "$yaml_file"
}

# Extract key paths from config
OUTPUT_DIR=$(parse_yaml_nested "$CONFIG_FILE" "paths" "output_dir")
LOG_DIR=$(parse_yaml_nested "$CONFIG_FILE" "paths" "log_dir")
SAMPLE_SHEET=$(parse_yaml_nested "$CONFIG_FILE" "paths" "sample_sheet")
PARTITION=$(parse_yaml_nested "$CONFIG_FILE" "slurm" "partition")
ACCOUNT=$(parse_yaml_nested "$CONFIG_FILE" "slurm" "account")

# Extract resource settings for each step
parse_resource() {
    local step="$1"
    local resource="$2"
    awk -v step="$step" -v res="$resource" '
        $0 ~ "^  "step":" { in_step=1; next }
        in_step && /^  [a-z]/ { in_step=0 }
        in_step && $0 ~ "^    "res":" {
            gsub(/.*: */, "");
            gsub(/["\047]/, "");
            gsub(/ *#.*$/, "");  # Remove comments
            print;
            exit
        }
    ' "$CONFIG_FILE"
}

# Preprocess resources
PREPROCESS_TIME=$(parse_resource "preprocess" "time")
PREPROCESS_MEM=$(parse_resource "preprocess" "mem")
PREPROCESS_CPUS=$(parse_resource "preprocess" "cpus")

# Demultiplex resources
DEMUX_TIME=$(parse_resource "demultiplex" "time")
DEMUX_MEM=$(parse_resource "demultiplex" "mem")
DEMUX_CPUS=$(parse_resource "demultiplex" "cpus")

# Quantify resources
QUANT_TIME=$(parse_resource "quantify" "time")
QUANT_MEM=$(parse_resource "quantify" "mem")
QUANT_CPUS=$(parse_resource "quantify" "cpus")

# Handle relative paths
if [[ ! "$SAMPLE_SHEET" = /* ]]; then
    SAMPLE_SHEET="${SCRIPT_DIR}/${SAMPLE_SHEET}"
fi

# Create output directories
FULL_LOG_DIR="${OUTPUT_DIR}/${LOG_DIR}"
mkdir -p "$OUTPUT_DIR" "$FULL_LOG_DIR"

echo "Output directory: $OUTPUT_DIR"
echo "Log directory: $FULL_LOG_DIR"
echo "Sample sheet: $SAMPLE_SHEET"
echo ""

# Validate sample sheet exists
if [[ ! -f "$SAMPLE_SHEET" ]]; then
    echo "ERROR: Sample sheet not found: $SAMPLE_SHEET"
    exit 1
fi

# Count samples
TOTAL_SAMPLES=$(tail -n +2 "$SAMPLE_SHEET" | grep -c "^0," || echo "0")
SKIPPED_SAMPLES=$(tail -n +2 "$SAMPLE_SHEET" | grep -c "^1," || echo "0")
echo "Total samples: $TOTAL_SAMPLES (skipped: $SKIPPED_SAMPLES)"
echo ""

# -----------------------------------------------------------------------------
# Build SBATCH common options
# -----------------------------------------------------------------------------
SBATCH_OPTS=""
if [[ -n "$PARTITION" && "$PARTITION" != "default" ]]; then
    SBATCH_OPTS="$SBATCH_OPTS --partition=$PARTITION"
fi
if [[ -n "$ACCOUNT" ]]; then
    SBATCH_OPTS="$SBATCH_OPTS --account=$ACCOUNT"
fi

# -----------------------------------------------------------------------------
# Submit jobs with dependencies
# -----------------------------------------------------------------------------
echo "Submitting jobs..."
echo ""
echo "Resource allocation from config:"
echo "  Preprocess:   ${PREPROCESS_CPUS} CPUs, ${PREPROCESS_MEM} memory, ${PREPROCESS_TIME} time"
echo "  Demultiplex:  ${DEMUX_CPUS} CPUs, ${DEMUX_MEM} memory, ${DEMUX_TIME} time"
echo "  Quantify:     ${QUANT_CPUS} CPUs, ${QUANT_MEM} memory, ${QUANT_TIME} time"
echo ""

# Job 1: Preprocess (trim + combine barcodes)
echo "Submitting 01_preprocess..."
JOB1=$(sbatch $SBATCH_OPTS \
    --cpus-per-task="${PREPROCESS_CPUS}" \
    --mem="${PREPROCESS_MEM}" \
    --time="${PREPROCESS_TIME}" \
    --output="${FULL_LOG_DIR}/01_preprocess_%j.out" \
    --error="${FULL_LOG_DIR}/01_preprocess_%j.err" \
    --export=CONFIG_FILE="$CONFIG_FILE",PIPELINE_DIR="$SCRIPT_DIR" \
    "${SCRIPT_DIR}/scripts/01_preprocess.slurm" | awk '{print $4}')
echo "  Job ID: $JOB1"

# Job 2: Demultiplex (depends on preprocess)
echo "Submitting 02_demultiplex..."
JOB2=$(sbatch $SBATCH_OPTS \
    --dependency=afterok:$JOB1 \
    --cpus-per-task="${DEMUX_CPUS}" \
    --mem="${DEMUX_MEM}" \
    --time="${DEMUX_TIME}" \
    --output="${FULL_LOG_DIR}/02_demultiplex_%j.out" \
    --error="${FULL_LOG_DIR}/02_demultiplex_%j.err" \
    --export=CONFIG_FILE="$CONFIG_FILE",PIPELINE_DIR="$SCRIPT_DIR" \
    "${SCRIPT_DIR}/scripts/02_demultiplex.slurm" | awk '{print $4}')
echo "  Job ID: $JOB2 (depends on $JOB1)"

# Job 3: Quantify (depends on demultiplex)
echo "Submitting 03_quantify..."
JOB3=$(sbatch $SBATCH_OPTS \
    --dependency=afterok:$JOB2 \
    --cpus-per-task="${QUANT_CPUS}" \
    --mem="${QUANT_MEM}" \
    --time="${QUANT_TIME}" \
    --output="${FULL_LOG_DIR}/03_quantify_%j.out" \
    --error="${FULL_LOG_DIR}/03_quantify_%j.err" \
    --export=CONFIG_FILE="$CONFIG_FILE",PIPELINE_DIR="$SCRIPT_DIR" \
    "${SCRIPT_DIR}/scripts/03_quantify.slurm" | awk '{print $4}')
echo "  Job ID: $JOB3 (depends on $JOB2)"

echo ""
echo "=============================================="
echo "Pipeline submitted successfully!"
echo "=============================================="
echo ""
echo "Job chain: $JOB1 -> $JOB2 -> $JOB3"
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER"
echo ""
echo "Check logs in: $FULL_LOG_DIR"
echo ""
