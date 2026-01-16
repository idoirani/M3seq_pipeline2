#!/bin/bash
# =============================================================================
# Debug script for plate barcode demultiplexing
# =============================================================================
# Usage: ./debug_plate_barcodes.sh <plate_bc_file> <sample_sheet>
#
# This script helps diagnose why plate barcode demultiplexing isn't matching.
# =============================================================================

set -eo pipefail

PLATE_BC_FILE="${1:-/scratch/gpfs/GITAI/ido/M3seq/output_scifi09_pipe2/02_demultiplexed/i7_demux/plate_bc/i7_AAGAGGCA.fastq.gz}"
SAMPLE_SHEET="${2:-/Users/ido/Dropbox/Bio/M3seq/Ziang_pipeline/pipeline/samples.csv}"

echo "=============================================="
echo "Plate Barcode Debug Analysis"
echo "=============================================="
echo ""
echo "Input file: $PLATE_BC_FILE"
echo "Sample sheet: $SAMPLE_SHEET"
echo ""

# Check if file exists
if [[ ! -f "$PLATE_BC_FILE" ]]; then
    echo "ERROR: Plate BC file not found: $PLATE_BC_FILE"
    exit 1
fi

# Show first 20 reads from the plate_bc file
echo "=============================================="
echo "First 10 reads from plate_bc file:"
echo "=============================================="
zcat "$PLATE_BC_FILE" 2>/dev/null | head -40 || gzip -dc "$PLATE_BC_FILE" | head -40

echo ""
echo "=============================================="
echo "First 13bp (plate barcode region) of first 100 reads:"
echo "=============================================="
zcat "$PLATE_BC_FILE" 2>/dev/null | awk 'NR%4==2 {print substr($0, 1, 13)}' | head -100 | sort | uniq -c | sort -rn | head -20

echo ""
echo "=============================================="
echo "Expected plate barcodes from sample sheet (first 20):"
echo "=============================================="
tail -n +2 "$SAMPLE_SHEET" | cut -d',' -f4 | sort | uniq | head -20

echo ""
echo "=============================================="
echo "Reverse complement of expected plate barcodes (first 20):"
echo "=============================================="
tail -n +2 "$SAMPLE_SHEET" | cut -d',' -f4 | sort | uniq | head -20 | while read bc; do
    echo "$bc" | tr 'ACGT' 'TGCA' | rev
done

echo ""
echo "=============================================="
echo "Read structure check - first 5 full reads:"
echo "=============================================="
echo ""
echo "Read structure (after UMI trimming):"
echo "  Position 1-13: Plate barcode"
echo "  Position 14-18: Skip region"
echo ""

zcat "$PLATE_BC_FILE" 2>/dev/null | head -20 | awk '
NR%4==1 {name=$0}
NR%4==2 {
    seq=$0
    print name
    printf "  Full:   %s (%d bp)\n", seq, length(seq)
    printf "  BC:     %s (pos 1-13)\n", substr(seq, 1, 13)
    printf "  Skip:   %s (pos 14-18)\n", substr(seq, 14, 5)
    print ""
}'

echo "=============================================="
echo "Barcode match summary"
echo "=============================================="
echo ""
echo "Checking if any of the expected barcodes appear in the first 10000 reads..."

# Get expected barcodes
EXPECTED=$(tail -n +2 "$SAMPLE_SHEET" | cut -d',' -f4 | sort | uniq)

# Get actual barcodes (first 13bp) from first 10000 reads
ACTUAL=$(zcat "$PLATE_BC_FILE" 2>/dev/null | awk 'NR%4==2 && NR<=40000 {print substr($0, 1, 13)}' | sort | uniq)

# Count matches
MATCH_COUNT=0
for bc in $EXPECTED; do
    if echo "$ACTUAL" | grep -q "^${bc}$"; then
        MATCH_COUNT=$((MATCH_COUNT + 1))
    fi
done

echo "Expected unique barcodes: $(echo "$EXPECTED" | wc -l)"
echo "Actual unique barcodes (first 10k reads): $(echo "$ACTUAL" | wc -l)"
echo "Matching barcodes: $MATCH_COUNT"

# Also check reverse complements
echo ""
echo "Checking reverse complement matches..."
RC_MATCH_COUNT=0
for bc in $EXPECTED; do
    rc=$(echo "$bc" | tr 'ACGT' 'TGCA' | rev)
    if echo "$ACTUAL" | grep -q "^${rc}$"; then
        RC_MATCH_COUNT=$((RC_MATCH_COUNT + 1))
    fi
done
echo "Matching reverse complement barcodes: $RC_MATCH_COUNT"

echo ""
echo "=============================================="
echo "Recommendation"
echo "=============================================="
if [[ $RC_MATCH_COUNT -gt $MATCH_COUNT ]]; then
    echo "*** Reverse complement barcodes match better! ***"
    echo "The plate barcodes in the sample sheet need to be reverse complemented."
    echo ""
    echo "To fix: Update the Python barcode generator in 02_demultiplex.slurm"
    echo "to reverse complement BARCODE_1 before writing to the FASTA file."
elif [[ $MATCH_COUNT -eq 0 && $RC_MATCH_COUNT -eq 0 ]]; then
    echo "*** Neither forward nor reverse complement barcodes match! ***"
    echo "Possible issues:"
    echo "  1. Wrong position - barcodes may not start at position 1"
    echo "  2. Different barcode length"
    echo "  3. Sequencing error or quality issues"
    echo "  4. Wrong sample sheet"
else
    echo "Forward barcodes match. The demultiplexing should be working."
    echo "Check cutadapt error rate settings or read quality."
fi
