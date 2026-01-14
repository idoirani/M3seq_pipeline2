#!/bin/bash
# =============================================================================
# 03_quantify.slurm - Quantification and statistics
# =============================================================================
# This script:
#   1. Optionally runs fastp QC on transcript reads
#   2. Runs quantification (kallisto-bustools or alevin-fry)
#   3. Generates per-library statistics
# =============================================================================

#SBATCH --job-name=m3seq_quantify
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00

# Note: Actual resource values are overridden by run_pipeline.sh from config.yaml

set -euo pipefail

# -----------------------------------------------------------------------------
# Load configuration
# -----------------------------------------------------------------------------
echo "=============================================="
echo "M3-seq Quantification"
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
elif isinstance(value, bool):
    print(str(value).lower())
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
conda activate "$CONDA_ENV"

# -----------------------------------------------------------------------------
# Get paths and settings from config
# -----------------------------------------------------------------------------
OUTPUT_DIR=$(parse_yaml "paths.output_dir")
WHITELIST=$(parse_yaml "paths.whitelist")
GENOME_DIR=$(parse_yaml "paths.genome_dir")

# QC settings
QC_ENABLED=$(parse_yaml "qc.enabled")

# Quantification settings
QUANT_TOOL=$(parse_yaml "quantification.tool")

MERGED_DIR="${OUTPUT_DIR}/03_merged_libraries"
QC_DIR="${OUTPUT_DIR}/04_qc"
QUANT_DIR="${OUTPUT_DIR}/05_quantification"
STATS_DIR="${OUTPUT_DIR}/06_stats"

mkdir -p "$QC_DIR" "$QUANT_DIR" "$STATS_DIR"

echo "Merged libraries directory: $MERGED_DIR"
echo "Quantification tool: $QUANT_TOOL"
echo "QC enabled: $QC_ENABLED"
echo ""

# -----------------------------------------------------------------------------
# Get list of libraries to process
# -----------------------------------------------------------------------------
LIBRARIES=$(ls "${MERGED_DIR}"/*.barcode_umi.fastq.gz 2>/dev/null | xargs -I{} basename {} .barcode_umi.fastq.gz)

if [[ -z "$LIBRARIES" ]]; then
    echo "ERROR: No merged library files found in $MERGED_DIR"
    exit 1
fi

echo "Libraries to process:"
for LIB in $LIBRARIES; do
    echo "  - $LIB"
done
echo ""

# -----------------------------------------------------------------------------
# Step 1: Optional QC with fastp
# -----------------------------------------------------------------------------
if [[ "$QC_ENABLED" == "true" ]]; then
    echo "=============================================="
    echo "Step 1: Running fastp QC"
    echo "=============================================="

    TRIM_POLY_G=$(parse_yaml "qc.fastp_options.trim_poly_g")
    QUAL_PHRED=$(parse_yaml "qc.fastp_options.qualified_quality_phred")
    LEN_REQ=$(parse_yaml "qc.fastp_options.length_required")

    for LIB in $LIBRARIES; do
        echo "Processing library: $LIB"

        BC_IN="${MERGED_DIR}/${LIB}.barcode_umi.fastq.gz"
        TX_IN="${MERGED_DIR}/${LIB}.transcript.fastq.gz"
        BC_OUT="${QC_DIR}/${LIB}.barcode_umi.qc.fastq.gz"
        TX_OUT="${QC_DIR}/${LIB}.transcript.qc.fastq.gz"

        FASTP_OPTS=""
        if [[ "$TRIM_POLY_G" == "true" ]]; then
            FASTP_OPTS="$FASTP_OPTS --trim_poly_g"
        fi

        # Run fastp on paired reads (barcode + transcript)
        fastp -i "$TX_IN" -I "$BC_IN" \
              -o "$TX_OUT" -O "$BC_OUT" \
              --qualified_quality_phred "$QUAL_PHRED" \
              --length_required "$LEN_REQ" \
              $FASTP_OPTS \
              --json "${QC_DIR}/${LIB}.fastp.json" \
              --html "${QC_DIR}/${LIB}.fastp.html" \
              --thread 4

        echo "  QC complete for $LIB"
    done
    echo ""

    # Update input directory for quantification
    INPUT_DIR="$QC_DIR"
    BC_SUFFIX=".barcode_umi.qc.fastq.gz"
    TX_SUFFIX=".transcript.qc.fastq.gz"
else
    echo "Step 1: QC skipped (disabled in config)"
    echo ""
    INPUT_DIR="$MERGED_DIR"
    BC_SUFFIX=".barcode_umi.fastq.gz"
    TX_SUFFIX=".transcript.fastq.gz"
fi

# -----------------------------------------------------------------------------
# Step 2: Quantification
# -----------------------------------------------------------------------------
echo "=============================================="
echo "Step 2: Running quantification ($QUANT_TOOL)"
echo "=============================================="

if [[ "$QUANT_TOOL" == "kallisto" ]]; then
    # Kallisto-bustools (kb) quantification
    KB_INDEX="${GENOME_DIR}/$(parse_yaml 'paths.kallisto_index')"
    KB_T2G="${GENOME_DIR}/$(parse_yaml 'paths.kallisto_t2g')"
    KB_THREADS=$(parse_yaml "quantification.kallisto.threads")
    KB_CHEMISTRY=$(parse_yaml "quantification.kallisto.chemistry")
    KB_STRAND=$(parse_yaml "quantification.kallisto.strand")

    echo "Kallisto index: $KB_INDEX"
    echo "T2G file: $KB_T2G"
    echo "Chemistry: $KB_CHEMISTRY"
    echo ""

    for LIB in $LIBRARIES; do
        echo "Quantifying library: $LIB"

        BC_FILE="${INPUT_DIR}/${LIB}${BC_SUFFIX}"
        TX_FILE="${INPUT_DIR}/${LIB}${TX_SUFFIX}"
        OUT_DIR="${QUANT_DIR}/${LIB}"

        mkdir -p "$OUT_DIR"

        kb count \
            -i "$KB_INDEX" \
            -g "$KB_T2G" \
            -x "$KB_CHEMISTRY" \
            -o "$OUT_DIR" \
            -t "$KB_THREADS" \
            -w "$WHITELIST" \
            --strand "$KB_STRAND" \
            --h5ad \
            --report \
            --overwrite \
            "$BC_FILE" "$TX_FILE"

        echo "  Quantification complete for $LIB"
        echo ""
    done

elif [[ "$QUANT_TOOL" == "alevin" ]]; then
    # Alevin-fry (simpleaf) quantification
    AF_INDEX="${GENOME_DIR}/$(parse_yaml 'paths.alevin_index')"
    AF_T2G="${GENOME_DIR}/$(parse_yaml 'paths.alevin_t2g')"
    AF_THREADS=$(parse_yaml "quantification.alevin.threads")
    AF_CHEMISTRY=$(parse_yaml "quantification.alevin.chemistry")
    AF_RESOLUTION=$(parse_yaml "quantification.alevin.resolution")
    AF_ORI=$(parse_yaml "quantification.alevin.expected_ori")

    echo "Alevin index: $AF_INDEX"
    echo "T2G file: $AF_T2G"
    echo "Chemistry: $AF_CHEMISTRY"
    echo ""

    for LIB in $LIBRARIES; do
        echo "Quantifying library: $LIB"

        BC_FILE="${INPUT_DIR}/${LIB}${BC_SUFFIX}"
        TX_FILE="${INPUT_DIR}/${LIB}${TX_SUFFIX}"
        OUT_DIR="${QUANT_DIR}/${LIB}"

        mkdir -p "$OUT_DIR"

        simpleaf quant \
            --reads1 "$BC_FILE" \
            --reads2 "$TX_FILE" \
            --threads "$AF_THREADS" \
            --index "$AF_INDEX" \
            --chemistry "$AF_CHEMISTRY" \
            --resolution "$AF_RESOLUTION" \
            --expected-ori "$AF_ORI" \
            --t2g-map "$AF_T2G" \
            --unfiltered-pl "$WHITELIST" \
            --output "$OUT_DIR" \
            --anndata-out

        echo "  Quantification complete for $LIB"
        echo ""
    done

else
    echo "ERROR: Unknown quantification tool: $QUANT_TOOL"
    echo "Valid options: kallisto, alevin"
    exit 1
fi

# -----------------------------------------------------------------------------
# Step 3: Generate statistics
# -----------------------------------------------------------------------------
echo "=============================================="
echo "Step 3: Generating statistics"
echo "=============================================="

# Create summary statistics
STATS_FILE="${STATS_DIR}/pipeline_stats.tsv"
echo -e "library\ttotal_reads\tmapped_reads\tcells\tgenes\tmedian_umi\tmedian_genes" > "$STATS_FILE"

for LIB in $LIBRARIES; do
    echo "Collecting stats for: $LIB"

    # Get read counts from input
    BC_FILE="${INPUT_DIR}/${LIB}${BC_SUFFIX}"
    TOTAL_READS=$(zcat "$BC_FILE" | wc -l)
    TOTAL_READS=$((TOTAL_READS / 4))

    # Get stats from quantification output
    QUANT_OUT="${QUANT_DIR}/${LIB}"

    if [[ "$QUANT_TOOL" == "kallisto" ]]; then
        # Parse kallisto output
        if [[ -f "${QUANT_OUT}/run_info.json" ]]; then
            MAPPED_READS=$(python3 -c "import json; print(json.load(open('${QUANT_OUT}/run_info.json'))['n_pseudoaligned'])" 2>/dev/null || echo "NA")
        else
            MAPPED_READS="NA"
        fi

        # Count cells and genes from matrix
        if [[ -f "${QUANT_OUT}/counts_unfiltered/cells_x_genes.mtx" ]]; then
            MATRIX_INFO=$(head -3 "${QUANT_OUT}/counts_unfiltered/cells_x_genes.mtx" | tail -1)
            CELLS=$(echo "$MATRIX_INFO" | awk '{print $1}')
            GENES=$(echo "$MATRIX_INFO" | awk '{print $2}')
        else
            CELLS="NA"
            GENES="NA"
        fi

    elif [[ "$QUANT_TOOL" == "alevin" ]]; then
        # Parse alevin output
        MAPPED_READS="NA"  # Would need to parse alevin logs

        # Count from output
        if [[ -d "${QUANT_OUT}/af_quant" ]]; then
            CELLS=$(find "${QUANT_OUT}" -name "*.h5ad" -exec python3 -c "import scanpy as sc; adata=sc.read_h5ad('{}'); print(adata.n_obs)" \; 2>/dev/null || echo "NA")
            GENES=$(find "${QUANT_OUT}" -name "*.h5ad" -exec python3 -c "import scanpy as sc; adata=sc.read_h5ad('{}'); print(adata.n_vars)" \; 2>/dev/null || echo "NA")
        else
            CELLS="NA"
            GENES="NA"
        fi
    fi

    # Placeholder for median stats (would need full matrix parsing)
    MEDIAN_UMI="NA"
    MEDIAN_GENES="NA"

    echo -e "${LIB}\t${TOTAL_READS}\t${MAPPED_READS}\t${CELLS}\t${GENES}\t${MEDIAN_UMI}\t${MEDIAN_GENES}" >> "$STATS_FILE"
done

echo ""
echo "Statistics saved to: $STATS_FILE"
echo ""

# Print summary table
echo "=============================================="
echo "Pipeline Summary"
echo "=============================================="
column -t -s $'\t' "$STATS_FILE"
echo ""

# -----------------------------------------------------------------------------
# Generate detailed per-library reports
# -----------------------------------------------------------------------------
echo "Generating detailed reports..."

python3 << 'PYTHON_SCRIPT'
import os
import json
import csv

stats_dir = os.environ.get('STATS_DIR', '')
quant_dir = os.environ.get('QUANT_DIR', '')
quant_tool = os.environ.get('QUANT_TOOL', 'kallisto')

# Read the stats file
stats_file = os.path.join(stats_dir, 'pipeline_stats.tsv')
with open(stats_file, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    libraries = list(reader)

# Generate a summary report
report_file = os.path.join(stats_dir, 'pipeline_report.txt')
with open(report_file, 'w') as f:
    f.write("=" * 60 + "\n")
    f.write("M3-seq Pipeline Report\n")
    f.write("=" * 60 + "\n\n")

    total_reads = 0
    total_cells = 0

    for lib in libraries:
        f.write(f"Library: {lib['library']}\n")
        f.write("-" * 40 + "\n")
        f.write(f"  Total reads:   {lib['total_reads']:>12}\n")
        f.write(f"  Mapped reads:  {lib['mapped_reads']:>12}\n")
        f.write(f"  Cells:         {lib['cells']:>12}\n")
        f.write(f"  Genes:         {lib['genes']:>12}\n")
        f.write("\n")

        try:
            total_reads += int(lib['total_reads'])
        except:
            pass
        try:
            total_cells += int(lib['cells'])
        except:
            pass

    f.write("=" * 60 + "\n")
    f.write("Totals\n")
    f.write("=" * 60 + "\n")
    f.write(f"  Total reads across all libraries: {total_reads:,}\n")
    f.write(f"  Total cells across all libraries: {total_cells:,}\n")
    f.write("\n")

print(f"Report saved to: {report_file}")
PYTHON_SCRIPT

export STATS_DIR QUANT_DIR QUANT_TOOL

echo ""
echo "End time: $(date)"
echo "Quantification complete!"
