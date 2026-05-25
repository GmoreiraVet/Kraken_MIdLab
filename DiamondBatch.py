#!/usr/bin/env python3

import os
import glob
import subprocess

# =========================================================
# USER SETTINGS
# =========================================================

# Folder containing FASTQ/FASTQ.GZ files
INPUT_DIR = "/home/viroicbas2023/Documents/Gmoreira/Carraças_Metagenomica_Gui/Reads_DeRNAfied/Reads_No_rRna"

# Output directory for DIAMOND results
OUTPUT_DIR = "/home/viroicbas2023/Documents/Gmoreira/Carraças_Metagenomica_Gui/Reads_DeRNAfied/Diamond_output"

# DIAMOND database (.dmnd)
DIAMOND_DB = "/home/viroicbas2023/Documents/DATABASES/NCBI_VIRAL/Viral_NCBI_DIAMOND.dmnd"

# Number of threads
THREADS = 30

# File extension to search for
PATTERN = "*.fastq.gz"

# =========================================================
# FILTER SETTINGS
# =========================================================

MIN_IDENTITY   = 70.0   # Minimum % identity (0–100). Set to None to disable.
MAX_EVALUE     = 1e-5   # Maximum e-value. Set to None to disable.
MIN_BITSCORE   = 0     # Minimum bit score. Set to None to disable.
MIN_ALN_LEN    = 30     # Minimum alignment length (aa). Set to None to disable.
MIN_QUERY_COV  = 30.0   # Minimum query coverage % (0–100). Set to None to disable.
MIN_SUBJECT_COV = None  # Minimum subject coverage % (0–100). Set to None to disable.

# =========================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Find all FASTQ files
files = glob.glob(os.path.join(INPUT_DIR, PATTERN))

if len(files) == 0:
    print("No input files found.")
    exit()

for f in files:

    # Sample name
    sample = os.path.basename(f)

    # Remove extensions
    for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        sample = sample.replace(ext, "")

    # Output file
    output_file = os.path.join(
        OUTPUT_DIR,
        f"{sample}_diamond.tsv"
    )

    print(f"\nProcessing: {sample}")

    # DIAMOND command
    cmd = [
        "diamond",
        "blastx",
        "-d", DIAMOND_DB,
        "-q", f,
        "-o", output_file,
        "-f", "6",
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "evalue",
        "bitscore",
        "stitle",
        "--top", "5",
        "--sensitive",
        "-p", str(THREADS)
    ]

    # --- Filters (only added if not None) ---
    if MIN_IDENTITY is not None:
        cmd += ["--id", str(MIN_IDENTITY)]

    if MAX_EVALUE is not None:
        cmd += ["-e", str(MAX_EVALUE)]

    if MIN_BITSCORE is not None:
        cmd += ["--min-score", str(MIN_BITSCORE)]

    if MIN_ALN_LEN is not None:
        cmd += ["--min-orf", str(MIN_ALN_LEN)]   # min ORF length in aa for blastx

    if MIN_QUERY_COV is not None:
        cmd += ["--query-cover", str(MIN_QUERY_COV)]

    if MIN_SUBJECT_COV is not None:
        cmd += ["--subject-cover", str(MIN_SUBJECT_COV)]

    print("Running command:")
    print(" ".join(cmd))

    # Run DIAMOND
    result = subprocess.run(cmd)

    if result.returncode == 0:
        print(f"Finished: {sample}")
    else:
        print(f"ERROR processing: {sample}")

print("\nAll samples completed.")