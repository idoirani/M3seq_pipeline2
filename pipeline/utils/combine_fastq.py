#!/usr/bin/env python3
"""
combine_fastq.py - Combine two FASTQ files by concatenating sequences

This script takes two FASTQ files and combines them read-by-read, concatenating
the sequences and quality scores. Useful for combining cell barcode and UMI
reads into a single file for downstream processing.

Usage:
    python combine_fastq.py read1.fastq.gz read2.fastq.gz output.fastq.gz

The output will have:
    - Sequence: read1_seq + read2_seq
    - Quality:  read1_qual + read2_qual
    - Header:   preserved from read1
"""

import sys
import gzip
import argparse
from pathlib import Path


def open_fastq(path, mode='rt'):
    """Open a FASTQ file, handling gzip compression automatically."""
    if str(path).endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_fastq(handle):
    """Generator that yields FASTQ records as tuples (header, seq, qual)."""
    while True:
        header = handle.readline().strip()
        if not header:
            break
        seq = handle.readline().strip()
        plus = handle.readline().strip()
        qual = handle.readline().strip()
        yield (header, seq, qual)


def combine_fastq(file1, file2, output, buffer_size=100000):
    """
    Combine two FASTQ files by concatenating sequences.

    Args:
        file1: Path to first FASTQ file (e.g., cell barcode)
        file2: Path to second FASTQ file (e.g., UMI)
        output: Path to output FASTQ file
        buffer_size: Number of records to buffer before writing
    """
    records_processed = 0
    buffer = []

    with open_fastq(file1, 'rt') as f1, open_fastq(file2, 'rt') as f2:
        iter1 = parse_fastq(f1)
        iter2 = parse_fastq(f2)

        # Open output file
        if str(output).endswith('.gz'):
            out_handle = gzip.open(output, 'wt', compresslevel=6)
        else:
            out_handle = open(output, 'wt')

        try:
            for (h1, s1, q1), (h2, s2, q2) in zip(iter1, iter2):
                # Validate read IDs match (strip trailing /1 /2 or space-separated info)
                id1 = h1.split()[0].replace('/1', '').replace('/2', '')
                id2 = h2.split()[0].replace('/1', '').replace('/2', '')

                if id1 != id2:
                    raise ValueError(
                        f"Read ID mismatch at record {records_processed + 1}:\n"
                        f"  File 1: {h1}\n"
                        f"  File 2: {h2}"
                    )

                # Combine sequences and qualities
                combined_seq = s1 + s2
                combined_qual = q1 + q2

                # Buffer the record
                buffer.append(f"{h1}\n{combined_seq}\n+\n{combined_qual}\n")
                records_processed += 1

                # Write buffer when full
                if len(buffer) >= buffer_size:
                    out_handle.writelines(buffer)
                    buffer.clear()

                    if records_processed % 1000000 == 0:
                        print(f"  Processed {records_processed:,} reads...", file=sys.stderr)

            # Write remaining records
            if buffer:
                out_handle.writelines(buffer)

        finally:
            out_handle.close()

    return records_processed


def main():
    parser = argparse.ArgumentParser(
        description='Combine two FASTQ files by concatenating sequences',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Combine cell barcode and UMI reads
    python combine_fastq.py cellbc.fastq.gz umi.fastq.gz combined.fastq.gz

    # The output will have concatenated sequences:
    # cellbc_seq + umi_seq (e.g., 16bp + 8bp = 24bp)
        """
    )

    parser.add_argument('file1', type=Path,
                        help='First FASTQ file (e.g., cell barcode)')
    parser.add_argument('file2', type=Path,
                        help='Second FASTQ file (e.g., UMI)')
    parser.add_argument('output', type=Path,
                        help='Output combined FASTQ file')
    parser.add_argument('--buffer-size', type=int, default=100000,
                        help='Number of records to buffer (default: 100000)')

    args = parser.parse_args()

    # Validate input files exist
    for f in [args.file1, args.file2]:
        if not f.exists():
            print(f"ERROR: Input file not found: {f}", file=sys.stderr)
            sys.exit(1)

    # Create output directory if needed
    args.output.parent.mkdir(parents=True, exist_ok=True)

    print(f"Combining FASTQ files:", file=sys.stderr)
    print(f"  File 1: {args.file1}", file=sys.stderr)
    print(f"  File 2: {args.file2}", file=sys.stderr)
    print(f"  Output: {args.output}", file=sys.stderr)

    try:
        n_reads = combine_fastq(
            args.file1,
            args.file2,
            args.output,
            buffer_size=args.buffer_size
        )
        print(f"Successfully combined {n_reads:,} reads", file=sys.stderr)

    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
