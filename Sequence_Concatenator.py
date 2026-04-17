#!/usr/bin/env python3

import os
from Bio import SeqIO
import gzip

# =========================
# CONFIGURATION SECTION
# =========================

CONFIG = {
    # Input directory (recursive)
    "input_dir": "/home/viroicbas2023/Documents/Gmoreira/EnteroVirusCarraçasUropigeaCustom/Enterovirus/ABRIL_2026/Seqs",

    # Output file
    "output_file": "/home/viroicbas2023/Documents/Gmoreira/EnteroVirusCarraçasUropigeaCustom/Enterovirus/ABRIL_2026/Seqs/concatenated_sequences.fasta",

    # Force format: "fasta", "fastq", or None (auto-detect)
    "format": None,

    # File extensions to include
    "extensions": (".fasta", ".fa", ".fna", ".fastq", ".fq", ".gz"),

    # Recursive search
    "recursive": True,
}

# =========================
# HELPERS
# =========================

def open_file(filepath):
    """Handle gzipped or plain text files."""
    if filepath.endswith(".gz"):
        return gzip.open(filepath, "rt")
    return open(filepath, "r")


def detect_format(filepath):
    """Detect format from first character."""
    with open_file(filepath) as f:
        first_char = f.read(1)
        if first_char == ">":
            return "fasta"
        elif first_char == "@":
            return "fastq"
        else:
            raise ValueError(f"Cannot detect format for {filepath}")


def collect_files(input_dir, extensions, recursive=True):
    """Collect files from directory."""
    collected = []
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.lower().endswith(extensions):
                collected.append(os.path.join(root, file))
        if not recursive:
            break
    return sorted(collected)


# =========================
# MAIN LOGIC
# =========================

def main():
    files = collect_files(
        CONFIG["input_dir"],
        CONFIG["extensions"],
        CONFIG["recursive"]
    )

    if not files:
        raise RuntimeError("No valid sequence files found.")

    print(f"Found {len(files)} files.")

    # Determine format
    if CONFIG["format"]:
        fmt = CONFIG["format"]
    else:
        fmt = detect_format(files[0])

    print(f"Using format: {fmt}")

    total_records = 0

    with open(CONFIG["output_file"], "w") as out_handle:

        for file in files:
            print(f"Processing: {file}")

            file_fmt = detect_format(file)
            if file_fmt != fmt:
                raise ValueError(f"Mixed formats detected: {file} is {file_fmt}, expected {fmt}")

            with open_file(file) as handle:
                records = SeqIO.parse(handle, fmt)

                count = SeqIO.write(records, out_handle, fmt)
                total_records += count

    print(f"\nDone.")
    print(f"Total records written: {total_records}")
    print(f"Output file: {CONFIG['output_file']}")


if __name__ == "__main__":
    main()