import sys
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def main(trim3_path, trim1_path, output_path, buffer_size=100000):
    with gzip.open(trim3_path, "rt") as f3, gzip.open(trim1_path, "rt") as f1:
        r3_iter = SeqIO.parse(f3, "fastq")
        r1_iter = SeqIO.parse(f1, "fastq")

        buffer = []
        with open("temp_combined.fastq", "wt") as out:
            for idx, (r3, r1) in enumerate(zip(r3_iter, r1_iter), 1):
                id3 = r3.id.split()[0].replace('/1', '').replace('/2', '')
                id1 = r1.id.split()[0].replace('/1', '').replace('/2', '')
                if id3 != id1:
                    raise ValueError(f"Read ID mismatch at record {idx}: {r3.id} vs {r1.id}")

                combined = SeqRecord(
                    seq=r3.seq + r1.seq,
                    id=r3.id,
                    description=r3.description,
                    letter_annotations={
                        "phred_quality": r3.letter_annotations["phred_quality"] + r1.letter_annotations["phred_quality"]
                    }
                )
                buffer.append(combined)

                if len(buffer) >= buffer_size:
                    SeqIO.write(buffer, out, "fastq")
                    buffer.clear()

            if buffer:
                SeqIO.write(buffer, out, "fastq")

    # Compress in parallel
    import subprocess
    subprocess.run(["pigz", "-p", "4", "-f", "temp_combined.fastq"])
    subprocess.run(["mv", "temp_combined.fastq.gz", output_path])

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python combine_fastq_buffered.py R3.fastq.gz R1.fastq.gz output.fastq.gz")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3])
