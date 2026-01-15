#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
from typing import TextIO

# ===== EDIT THESE PATHS =====
FASTQ_PATH = "/home/viroicbas2023/Documents/Gmoreira/Serginho/RatHevFalahado/4GN65X_fastq/4GN65X_1_S76_Liver.fastq"
FASTA_PATH = "/home/viroicbas2023/Documents/Gmoreira/Serginho/RatHevFalahado/Fasta.fasta"
# ===========================

def open_maybe_gzip(path: str, mode: str = "rt") -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8")

def fastq_to_fasta(fastq_path: str, fasta_path: str) -> None:
    with open_maybe_gzip(fastq_path, "rt") as fin, open(fasta_path, "w", encoding="utf-8") as fout:
        line_num = 0

        for line in fin:
            line_num += 1
            line = line.rstrip("\n")
            mod = line_num % 4

            if mod == 1:  # Header
                if not line.startswith("@"):
                    raise ValueError(f"Invalid FASTQ format at line {line_num}")
                header = line[1:].split()[0]
                fout.write(f">{header}\n")

            elif mod == 2:  # Sequence
                fout.write(line + "\n")

        if line_num % 4 != 0:
            raise ValueError("FASTQ file is truncated or malformed")

    print("Conversion complete.")
    print(f"Input : {fastq_path}")
    print(f"Output: {fasta_path}")

if __name__ == "__main__":
    fastq_to_fasta(FASTQ_PATH, FASTA_PATH)
