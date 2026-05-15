#!/usr/bin/env python3
"""
run_unicycler.py — Unicycler assembly wrapper
Edit the PARAMETERS section below, then run:  python run_unicycler.py
"""

import subprocess
import sys
import os
from datetime import datetime

# =============================================================================
#  PARAMETERS — edit these
# =============================================================================

# --- Input reads -------------------------------------------------------------
SHORT_R1    = ""          # Illumina R1  (leave empty string if not used)
SHORT_R2    = ""          # Illumina R2  (leave empty string if not used)
LONG_READS  = "/home/viroicbas2023/Documents/Gmoreira/GuilhermeMoreira2026/Tacheng/May2026/ReadsTogether/tacheng_joined.fastq"

# --- Output ------------------------------------------------------------------
OUT_DIR = "/home/viroicbas2023/Documents/Gmoreira/GuilhermeMoreira2026/Tacheng/May2026/UniCyler"

# --- Assembly mode -----------------------------------------------------------
# Options: "conservative" | "normal" | "bold"
#   conservative  → cautious bridging, fewer errors, may leave gaps
#   normal        → balanced (recommended default)
#   bold          → aggressive bridging, more complete but riskier
MODE = "normal"

# --- Performance -------------------------------------------------------------
THREADS = 30

# --- Verbosity ---------------------------------------------------------------
# Options: 0 (silent) | 1 (moderate) | 2 (full/debug)
VERBOSITY = 1

# --- Min read length (long reads) -------------------------------------------
# Reads shorter than this are discarded before assembly.
# Increase if assembly quality is poor (e.g. 1000, 2000).
MIN_LONG_READ_LENGTH = 200

# --- Minimum component size --------------------------------------------------
# Discard assembly graph components smaller than this (bp).
MIN_COMPONENT_SIZE = 1000

# --- Minimum dead end size ---------------------------------------------------
# Remove dead-end sequences shorter than this (bp).
MIN_DEAD_END_SIZE = 1000

# --- Keep intermediate files -------------------------------------------------
# 0 = remove all,  1 = keep some,  2 = keep everything (useful for debugging)
KEEP = 0

# =============================================================================
#  (no need to edit below this line)
# =============================================================================

def main():
    print()
    print("=" * 50)
    print(f"  Unicycler assembly — {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 50)

    # Validate short reads: must be both or neither
    if bool(SHORT_R1) != bool(SHORT_R2):
        sys.exit("ERROR: Provide both SHORT_R1 and SHORT_R2, or leave both empty.")

    # Validate that at least one input is given
    if not SHORT_R1 and not LONG_READS:
        sys.exit("ERROR: No input reads provided. Set SHORT_R1/R2 and/or LONG_READS.")

    # Validate input files exist
    for label, path in [("SHORT_R1", SHORT_R1), ("SHORT_R2", SHORT_R2), ("LONG_READS", LONG_READS)]:
        if path and not os.path.isfile(path):
            sys.exit(f"ERROR: {label} file not found:\n  {path}")

    # Build command
    cmd = ["unicycler"]

    if SHORT_R1 and SHORT_R2:
        cmd += ["-1", SHORT_R1, "-2", SHORT_R2]

    if LONG_READS:
        cmd += ["-l", LONG_READS]

    cmd += [
        "--out",                    OUT_DIR,
        "--mode",                   MODE,
        "--threads",                str(THREADS),
        "--verbosity",              str(VERBOSITY),
        "--min_long_read_length",   str(MIN_LONG_READ_LENGTH),
        "--min_component_size",     str(MIN_COMPONENT_SIZE),
        "--min_dead_end_size",      str(MIN_DEAD_END_SIZE),
        "--keep",                   str(KEEP),
    ]

    os.makedirs(OUT_DIR, exist_ok=True)

    print()
    print("Command:")
    print("  " + " ".join(cmd))
    print()

    result = subprocess.run(cmd)

    print()
    print("=" * 50)
    if result.returncode == 0:
        print(f"  Done — {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"  Results in: {OUT_DIR}")
    else:
        print(f"  Unicycler exited with code {result.returncode}")
    print("=" * 50)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()