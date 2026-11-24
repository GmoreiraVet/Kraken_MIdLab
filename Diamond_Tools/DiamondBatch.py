#!/usr/bin/env python3

import os
import glob
import subprocess

# =========================================================
# USER SETTINGS
# =========================================================

# Folder containing FASTQ/FASTQ.GZ files
INPUT_DIR = "/home/viroicbas2023/Documents/Gmoreira/VIRUS_Others/PATOS_INFLENZA_PARAMYXO/Reads"

# Output directory for DIAMOND results
OUTPUT_DIR = "/home/viroicbas2023/Documents/Gmoreira/VIRUS_Others/PATOS_INFLENZA_PARAMYXO/DiamondClassification"

# DIAMOND database (.dmnd)
DIAMOND_DB = "/home/viroicbas2023/Documents/DATABASES/NCBI_VIRAL/Viral_NCBI_DIAMOND.dmnd"

# Number of threads
THREADS = 30

# File extension to search for
PATTERN = "*.fastq.gz"

# =========================================================
# HIT REPORTING SETTINGS
# =========================================================
# Use ONE of the two options below (set the other to None).
#
# TOP_PERCENT: report all hits within X% of the best bit score per read
#              (e.g. 5 = hits scoring >= 95% of the best). Number of hits is uncapped.
# MAX_TARGET_SEQS: report at most N hits per read (1 = best hit only).
#
# If both are set, TOP_PERCENT takes precedence (DIAMOND ignores
# --max-target-seqs when --top is used), so a warning is printed.

TOP_PERCENT     = None      # e.g. 1, 5, 10, or None
MAX_TARGET_SEQS = 1  # e.g. 1, 5, 25, or None

# =========================================================
# FILTER SETTINGS
# =========================================================

MIN_IDENTITY    = 60.0   # Minimum % identity (0–100). Set to None to disable.
MAX_EVALUE      = 1e-5   # Maximum e-value. Set to None to disable.
MIN_BITSCORE    = None      # Minimum bit score. Set to None to disable.
MIN_ORF_LEN     = 30     # Minimum ORF length (aa) used in translation (--min-orf). None to disable.
MIN_ALN_LEN     = None   # Minimum alignment length (aa), applied after DIAMOND. None to disable.
MIN_QUERY_COV   = 30.0   # Minimum query coverage % (0–100). Set to None to disable.
MIN_SUBJECT_COV = None   # Minimum subject coverage % (0–100). Set to None to disable.

# =========================================================

# Validate hit settings
if TOP_PERCENT is not None and MAX_TARGET_SEQS is not None:
    print("WARNING: both TOP_PERCENT and MAX_TARGET_SEQS are set. "
          "DIAMOND will ignore MAX_TARGET_SEQS; using TOP_PERCENT only.")
    MAX_TARGET_SEQS = None

if TOP_PERCENT is None and MAX_TARGET_SEQS is None:
    print("NOTE: neither TOP_PERCENT nor MAX_TARGET_SEQS is set; "
          "DIAMOND default (--max-target-seqs 25) will be used.")

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

    # Output files
    output_file = os.path.join(OUTPUT_DIR, f"{sample}_diamond.tsv")
    tmp_file = output_file + ".tmp"

    print(f"\nProcessing: {sample}")

    # DIAMOND command
    cmd = [
        "diamond",
        "blastx",
        "-d", DIAMOND_DB,
        "-q", f,
        "-o", tmp_file if MIN_ALN_LEN is not None else output_file,
        "-f", "6",
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "evalue",
        "bitscore",
        "stitle",
        "--sensitive",
        "-p", str(THREADS)
    ]

    # --- Hit reporting ---
    if TOP_PERCENT is not None:
        cmd += ["--top", str(TOP_PERCENT)]
    elif MAX_TARGET_SEQS is not None:
        cmd += ["--max-target-seqs", str(MAX_TARGET_SEQS)]

    # --- Filters (only added if not None) ---
    if MIN_IDENTITY is not None:
        cmd += ["--id", str(MIN_IDENTITY)]

    if MAX_EVALUE is not None:
        cmd += ["-e", str(MAX_EVALUE)]

    if MIN_BITSCORE is not None:
        cmd += ["--min-score", str(MIN_BITSCORE)]

    if MIN_ORF_LEN is not None:
        cmd += ["--min-orf", str(MIN_ORF_LEN)]

    if MIN_QUERY_COV is not None:
        cmd += ["--query-cover", str(MIN_QUERY_COV)]

    if MIN_SUBJECT_COV is not None:
        cmd += ["--subject-cover", str(MIN_SUBJECT_COV)]

    print("Running command:")
    print(" ".join(cmd))

    # Run DIAMOND
    result = subprocess.run(cmd)

    if result.returncode != 0:
        print(f"ERROR processing: {sample}")
        continue

    # Post-filter on alignment length (column 4 = "length")
    if MIN_ALN_LEN is not None:
        kept = 0
        with open(tmp_file) as fin, open(output_file, "w") as fout:
            for line in fin:
                cols = line.rstrip("\n").split("\t")
                if len(cols) > 3 and int(cols[3]) >= MIN_ALN_LEN:
                    fout.write(line)
                    kept += 1
        os.remove(tmp_file)
        print(f"Alignment length filter (>= {MIN_ALN_LEN} aa): {kept} hits kept")

    print(f"Finished: {sample}")

print("\nAll samples completed.")