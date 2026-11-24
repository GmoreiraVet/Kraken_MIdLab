#!/usr/bin/env python3
"""
run_diamond_best_hit.py

Runs DIAMOND (blastx or blastp) restricted to the single best hit per
query read (--max-target-seqs 1), using the same tabular output format
expected by extract_species_reads.py.

No command-line arguments — edit the CONFIG block below and run:
    python run_diamond_best_hit.py
"""

import subprocess
import shutil
import sys

# ============================== CONFIG ==============================
DIAMOND_MODE = "blastx"          # "blastx" (nucleotide query) or "blastp" (protein query)

QUERY_FILE = "ReaDS/PL8XTX_7_S76_Liver.fastq"       # input reads (fasta/fastq, .fastq.gz also supported by DIAMOND)
DATABASE = "/home/viroicbas2023/Documents/DATABASES/DIAMOND/DIAMOND_RVDB.dmnd"  # path to the DIAMOND .dmnd database

OUTPUT_TSV = "/home/viroicbas2023/Documents/Andreia/ROCAHEPEVIRUS/rocas_DIAMOND.tsv"

EVALUE = "1e-10"
THREADS = 30

# outfmt 6 columns — must match what extract_species_reads.py expects
OUTFMT_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "qlen", "slen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle",
]
# ======================================================================


def main():
    diamond_bin = shutil.which("diamond")
    if diamond_bin is None:
        sys.exit("ERROR: 'diamond' executable not found in PATH. "
                  "Activate the correct conda env or check your installation.")

    cmd = [
        diamond_bin, DIAMOND_MODE,
        "--query", QUERY_FILE,
        "--db", DATABASE,
        "--out", OUTPUT_TSV,
        "--outfmt", "6", *OUTFMT_COLUMNS,
        "--max-target-seqs", "1",   # single best hit per read
        "--evalue", EVALUE,
        "--threads", str(THREADS),
    ]

    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd)

    if result.returncode != 0:
        sys.exit(f"DIAMOND exited with code {result.returncode}")

    print(f"Done. Best-hit-per-read results -> {OUTPUT_TSV}")


if __name__ == "__main__":
    main()