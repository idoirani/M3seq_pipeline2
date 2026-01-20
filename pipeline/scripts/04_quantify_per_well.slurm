#!/bin/bash
# =============================================================================
# 04_quantify_per_well.slurm - Quantification per well (preserving plate info)
# =============================================================================
# This script:
#   1. Runs kallisto quantification on EACH WELL separately (not merged)
#   2. Combines all outputs into a single cell index and expression matrix
# =============================================================================

#SBATCH --job-name=m3seq_quant_well
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00

set -eo pipefail

export PS1=${PS1:-""}

# -----------------------------------------------------------------------------
# Load configuration
# -----------------------------------------------------------------------------
echo "=============================================="
echo "M3-seq Per-Well Quantification"
echo "=============================================="
echo "Start time: $(date)"

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
source activate "$CONDA_ENV" || conda activate "$CONDA_ENV"

# Verify commands
command -v kb >/dev/null 2>&1 || { echo "ERROR: kb not found"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found"; exit 1; }

# -----------------------------------------------------------------------------
# Get paths and settings
# -----------------------------------------------------------------------------
OUTPUT_DIR=$(parse_yaml "paths.output_dir")
WHITELIST=$(parse_yaml "paths.whitelist")
GENOME_DIR=$(parse_yaml "paths.genome_dir")
SAMPLE_SHEET=$(parse_yaml "paths.sample_sheet")

# Handle relative sample sheet path
if [[ ! "$SAMPLE_SHEET" = /* ]]; then
    SAMPLE_SHEET="${PIPELINE_DIR}/${SAMPLE_SHEET}"
fi

KB_INDEX="${GENOME_DIR}/$(parse_yaml 'paths.kallisto_index')"
KB_T2G="${GENOME_DIR}/$(parse_yaml 'paths.kallisto_t2g')"
KB_THREADS=$(parse_yaml "quantification.kallisto.threads")
KB_CHEMISTRY=$(parse_yaml "quantification.kallisto.chemistry")
KB_STRAND=$(parse_yaml "quantification.kallisto.strand")

DEMUX_DIR="${OUTPUT_DIR}/02_demultiplexed"
PLATE_DEMUX_DIR="${DEMUX_DIR}/plate_demux"
WELL_QUANT_DIR="${OUTPUT_DIR}/07_well_quantification"
COMBINED_DIR="${OUTPUT_DIR}/08_combined_counts"

mkdir -p "$WELL_QUANT_DIR" "$COMBINED_DIR"

echo "Kallisto index: $KB_INDEX"
echo "T2G file: $KB_T2G"
echo "Output: $WELL_QUANT_DIR"
echo ""
echo "DEBUG: DEMUX_DIR = $DEMUX_DIR"
echo "DEBUG: PLATE_DEMUX_DIR = $PLATE_DEMUX_DIR"
echo "DEBUG: Sample sheet = $SAMPLE_SHEET"
echo ""

# -----------------------------------------------------------------------------
# Step 1: Run kallisto on each well
# -----------------------------------------------------------------------------
echo "=============================================="
echo "Step 1: Quantifying each well separately"
echo "=============================================="

# Get list of all samples from sample sheet (tr -d '\r' removes Windows carriage returns)
SAMPLES=$(tail -n +2 "$SAMPLE_SHEET" | grep "^0," | cut -d',' -f3 | tr -d '\r')
TOTAL_SAMPLES=$(echo "$SAMPLES" | wc -l)
CURRENT=0

for SAMPLE in $SAMPLES; do
    CURRENT=$((CURRENT + 1))
    echo ""
    echo "[$CURRENT/$TOTAL_SAMPLES] Processing: $SAMPLE"

    # Get i7 barcode for this sample (tr -d '\r' removes Windows carriage returns)
    I7_BC=$(grep ",${SAMPLE}," "$SAMPLE_SHEET" | cut -d',' -f5 | tr -d '\r')
    I7_NAME="i7_${I7_BC}"

    # Input files
    BC_FILE="${PLATE_DEMUX_DIR}/${I7_NAME}/barcode_umi/${SAMPLE}.fastq.gz"
    TX_FILE="${PLATE_DEMUX_DIR}/${I7_NAME}/transcript/${SAMPLE}.fastq.gz"

    echo "  DEBUG: I7_BC=$I7_BC, I7_NAME=$I7_NAME"
    echo "  DEBUG: BC_FILE=$BC_FILE"

    if [[ ! -f "$BC_FILE" ]]; then
        echo "  WARNING: Barcode file not found: $BC_FILE"
        continue
    fi
    if [[ ! -f "$TX_FILE" ]]; then
        echo "  WARNING: Transcript file not found: $TX_FILE"
        continue
    fi

    OUT_DIR="${WELL_QUANT_DIR}/${SAMPLE}"

    if [[ $CURRENT -le 3 ]]; then
        echo "  DEBUG bash: OUT_DIR=$OUT_DIR"
        echo "  DEBUG bash: SAMPLE='$SAMPLE'"
    fi

    # Skip if already processed
    if [[ -f "${OUT_DIR}/counts_unfiltered/cells_x_genes.mtx" && -f "${OUT_DIR}/raw_read_counts.tsv" ]]; then
        echo "  Already processed, skipping..."
        continue
    fi

    mkdir -p "$OUT_DIR"

    # Step 1a: Count raw reads per r2 barcode (before UMI deduplication)
    echo "  Counting raw reads per cell barcode..."
    zcat "$BC_FILE" | awk 'NR%4==2 {print substr($0, 1, 16)}' | sort | uniq -c | \
        awk '{print $2"\t"$1}' > "${OUT_DIR}/raw_read_counts.tsv"

    # Step 1b: Run kallisto
    echo "  Running kallisto..."
    kb count \
        -i "$KB_INDEX" \
        -g "$KB_T2G" \
        -x "$KB_CHEMISTRY" \
        -o "$OUT_DIR" \
        -t "$KB_THREADS" \
        -w "$WHITELIST" \
        --strand "$KB_STRAND" \
        --overwrite \
        "$BC_FILE" "$TX_FILE" \
        2>&1 | tail -5

    echo "  Done"
done

echo ""
echo "Per-well quantification complete."
echo ""

# -----------------------------------------------------------------------------
# Step 2: Combine all wells into single matrix
# -----------------------------------------------------------------------------
echo "=============================================="
echo "Step 2: Combining all wells into single matrix"
echo "=============================================="

export SAMPLE_SHEET WELL_QUANT_DIR COMBINED_DIR KB_T2G

python3 << 'PYTHON_SCRIPT'
import os
import gzip
import csv
import numpy as np
from scipy import sparse
from scipy.io import mmread
from collections import defaultdict

sample_sheet = os.environ['SAMPLE_SHEET']
well_quant_dir = os.environ['WELL_QUANT_DIR']
combined_dir = os.environ['COMBINED_DIR']
t2g_file = os.environ['KB_T2G']

print(f"DEBUG Python: well_quant_dir = '{well_quant_dir}'")
print(f"DEBUG Python: sample_sheet = '{sample_sheet}'")

print("Loading sample information...")

# Load sample sheet
samples = []
with open(sample_sheet, 'r', encoding='utf-8-sig') as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row['SKIP'] == '0':
            samples.append(row)

print(f"Found {len(samples)} samples")

# Extract well from sample name (last part after final underscore)
def get_well(sample_name):
    # e.g., "Lib1_EZ_Tn5_MG1655_Exponential_A1" -> "A1"
    return sample_name.split('_')[-1]

# Extract treatment (everything except the well)
def get_treatment(sample_name):
    # e.g., "Lib1_EZ_Tn5_MG1655_Exponential_A1" -> "Lib1_EZ_Tn5_MG1655_Exponential"
    parts = sample_name.split('_')
    return '_'.join(parts[:-1])

# Load gene list from t2g file
print("Loading gene list...")
genes = []
gene_to_idx = {}
with open(t2g_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            gene = parts[1]  # gene name is second column
            if gene not in gene_to_idx:
                gene_to_idx[gene] = len(genes)
                genes.append(gene)

print(f"Found {len(genes)} genes")

# Process each sample and collect data
print("Processing samples...")
all_cell_data = []
all_matrices = []
global_cell_idx = 0

for i, sample_info in enumerate(samples):
    sample_name = sample_info['SAMPLE_NAME']
    library_name = sample_info['LIBRARY_NAME']
    r1_barcode = sample_info['BARCODE_1']  # plate barcode
    i7_barcode = sample_info['BARCODE_2']  # i7 barcode

    well = get_well(sample_name)
    treatment = get_treatment(sample_name)

    # Paths to kallisto output
    mtx_path = os.path.join(well_quant_dir, sample_name, 'counts_unfiltered', 'cells_x_genes.mtx')
    bc_path = os.path.join(well_quant_dir, sample_name, 'counts_unfiltered', 'cells_x_genes.barcodes.txt')
    reads_path = os.path.join(well_quant_dir, sample_name, 'raw_read_counts.tsv')

    if not os.path.exists(mtx_path):
        if i < 3:  # Debug first 3 samples
            print(f"  DEBUG: Looking for mtx at: {mtx_path}")
            print(f"  DEBUG: sample_name = '{sample_name}'")
        print(f"  [{i+1}/{len(samples)}] {sample_name}: No output found, skipping")
        continue

    # Load matrix
    mtx = mmread(mtx_path).tocsr()

    # Load barcodes (r2 barcodes / 10x barcodes)
    with open(bc_path, 'r') as f:
        r2_barcodes = [line.strip() for line in f]

    # Load raw read counts per barcode
    raw_reads = {}
    if os.path.exists(reads_path):
        with open(reads_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    raw_reads[parts[0]] = int(parts[1])

    n_cells = mtx.shape[0]

    if n_cells == 0:
        print(f"  [{i+1}/{len(samples)}] {sample_name}: No cells, skipping")
        continue

    # Calculate per-cell stats
    umis_per_cell = np.array(mtx.sum(axis=1)).flatten()
    genes_per_cell = np.array((mtx > 0).sum(axis=1)).flatten()

    # Create cell index entries
    for j, r2_bc in enumerate(r2_barcodes):
        cell_id = f"{r2_bc}::{well}::{treatment}"
        cell_data = {
            'cell_index': global_cell_idx,
            'cell_id': cell_id,
            'r2_barcode': r2_bc,
            'r1_barcode': r1_barcode,
            'plate_well': well,
            'i7_barcode': i7_barcode,
            'treatment': treatment,
            'library_name': library_name,
            'sample_name': sample_name,
            'n_reads': raw_reads.get(r2_bc, 0),
            'n_umis': int(umis_per_cell[j]),
            'n_genes': int(genes_per_cell[j]),
        }
        all_cell_data.append(cell_data)
        global_cell_idx += 1

    all_matrices.append(mtx)
    print(f"  [{i+1}/{len(samples)}] {sample_name}: {n_cells} cells")

print(f"\nTotal cells: {global_cell_idx}")

# Combine all matrices
print("Combining matrices...")
if all_matrices:
    combined_matrix = sparse.vstack(all_matrices, format='csr')
    print(f"Combined matrix shape: {combined_matrix.shape}")
else:
    print("ERROR: No matrices to combine!")
    exit(1)

# Save gene list
gene_list_path = os.path.join(combined_dir, 'gene_list.txt')
print(f"Saving gene list to {gene_list_path}")
with open(gene_list_path, 'w') as f:
    for gene in genes:
        f.write(f"{gene}\n")

# Save cell index
cell_index_path = os.path.join(combined_dir, 'cell_index.tsv.gz')
print(f"Saving cell index to {cell_index_path}")
with gzip.open(cell_index_path, 'wt') as f:
    writer = csv.DictWriter(f, fieldnames=[
        'cell_index', 'cell_id', 'r2_barcode', 'r1_barcode',
        'plate_well', 'i7_barcode', 'treatment', 'library_name',
        'sample_name', 'n_reads', 'n_umis', 'n_genes'
    ], delimiter='\t')
    writer.writeheader()
    writer.writerows(all_cell_data)

# Save expression matrix as npz
matrix_path = os.path.join(combined_dir, 'expression_matrix.npz')
print(f"Saving expression matrix to {matrix_path}")
sparse.save_npz(matrix_path, combined_matrix)

# Print summary
print("\n" + "=" * 60)
print("Summary")
print("=" * 60)
print(f"Total cells: {len(all_cell_data):,}")
print(f"Total genes: {len(genes):,}")
print(f"Matrix shape: {combined_matrix.shape}")
print(f"Non-zero entries: {combined_matrix.nnz:,}")
print(f"Sparsity: {100 * (1 - combined_matrix.nnz / (combined_matrix.shape[0] * combined_matrix.shape[1])):.2f}%")
print("")
print("Output files:")
print(f"  {cell_index_path}")
print(f"  {matrix_path}")
print(f"  {gene_list_path}")

PYTHON_SCRIPT

echo ""
echo "End time: $(date)"
echo "Per-well quantification and combination complete!"
