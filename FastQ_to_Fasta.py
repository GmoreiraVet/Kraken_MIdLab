#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

# ===== EDIT THESE PATHS =====
INPUT_FOLDER = "/home/viroicbas2023/Downloads/YZ8NQ3_fastq"
OUTPUT_FOLDER = "/home/viroicbas2023/Downloads/YZ8NQ3_fastq/OutputFastaReruns"
TRIM_SIZE = 10
# ============================


def trim_sequence(seq: str, trim: int) -> str:
    """Trim bases from both ends of a sequence."""
    if len(seq) <= trim * 2:
        return ""
    return seq[trim:-trim]


def fastq_to_trimmed_fasta(input_file: Path, output_file: Path, trim: int):
    """Convert FASTQ to FASTA while trimming sequences."""

    with open(input_file, "r") as fin, open(output_file, "w") as fout:

        while True:
            header = fin.readline()
            if not header:
                break

            seq = fin.readline().strip()
            fin.readline()  # +
            fin.readline()  # quality

            header = header.strip()

            if not header.startswith("@"):
                raise ValueError(f"Invalid FASTQ header in {input_file}")

            trimmed = trim_sequence(seq, trim)

            if trimmed:
                read_id = header[1:].split()[0]
                fout.write(f">{read_id}\n{trimmed}\n")


def main():

    input_path = Path(INPUT_FOLDER)
    output_path = Path(OUTPUT_FOLDER)

    output_path.mkdir(parents=True, exist_ok=True)

    fastq_files = (
        list(input_path.glob("*.fastq")) +
        list(input_path.glob("*.fq"))
    )

    if not fastq_files:
        print("No FASTQ files found.")
        return

    for fastq in fastq_files:

        output_file = output_path / (fastq.stem + ".fasta")

        print(f"Processing: {fastq.name}")

        fastq_to_trimmed_fasta(fastq, output_file, TRIM_SIZE)

    print(f"\nFinished processing {len(fastq_files)} files.")
    print(f"Output folder: {output_path}")


if __name__ == "__main__":
    main()
