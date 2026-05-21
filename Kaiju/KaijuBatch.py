#!/usr/bin/env python3
"""
run_kaiju.py
------------
Runs Kaiju on all FASTA/FASTQ/GZ files in an input folder, then converts
each output file to a human-readable format with kaiju-addTaxonNames.

Edit the CONFIGURATION section below before running.
"""

import subprocess
import sys
from pathlib import Path

# =============================================================================
# CONFIGURATION — edit these paths and parameters before running
# =============================================================================

# --- Kaiju database files ---
FMI_FILE   = "/home/viroicbas2023/Documents/KaijuDB/kaiju_db_rvdb_2024-12-20/kaiju_db_rvdb.fmi"
NODES_FILE = "/home/viroicbas2023/Documents/KaijuDB/kaiju_db_rvdb_2024-12-20/nodes.dmp"
NAMES_FILE = "/home/viroicbas2023/Documents/KaijuDB/kaiju_db_rvdb_2024-12-20/names.dmp"

# --- Input / output folders ---
INPUT_FOLDER  = "/home/viroicbas2023/Documents/Gmoreira/AndreiaGuiMay8/KBFTRN_fastq"
OUTPUT_FOLDER = "/home/viroicbas2023/Documents/Gmoreira/AndreiaGuiMay8/KaijuResults_BIGDB"

# --- kaiju run parameters ---
KAIJU_THREADS    = 30       # -z  number of parallel threads
KAIJU_MODE       = "greedy" # -a  alignment mode: "mem" or "greedy"
KAIJU_MIN_SCORE  = 75       # -s  minimum match score (greedy mode)
KAIJU_MISMATCHES = 3        # -e  max mismatches (greedy mode)


# --- kaiju-addTaxonNames parameters ---
OMIT_UNCLASSIFIED_NAMES = True  # -u  omit unclassified reads from named output
PRINT_FULL_TAXON_PATH   = False  # -p  print full taxon path instead of just taxon name
TAXON_RANKS             = ""     # -r  comma-separated ranks to include in path, e.g. "phylum,genus"
                                 #     (only used when PRINT_FULL_TAXON_PATH is False)

# --- Paths to executables (leave as-is if they are on $PATH) ---
KAIJU_BIN           = "kaiju"
ADD_TAXON_NAMES_BIN = "kaiju-addTaxonNames"

# =============================================================================
# END OF CONFIGURATION
# =============================================================================

# Supported input extensions (checked from longest to shortest to handle .fastq.gz etc.)
SUPPORTED_EXTENSIONS = [
    ".fastq.gz", ".fasta.gz", ".fq.gz", ".fa.gz", ".fna.gz",
    ".fastq", ".fasta", ".fq", ".fa", ".fna",
]


def find_input_files(folder: str) -> list[Path]:
    """Return all FASTA/FASTQ/GZ files inside *folder* (non-recursive)."""
    p = Path(folder)
    if not p.is_dir():
        sys.exit(f"[ERROR] Input folder not found: {folder}")
    files = []
    for f in sorted(p.iterdir()):
        if f.is_file() and any(f.name.endswith(ext) for ext in SUPPORTED_EXTENSIONS):
            files.append(f)
    return files


def sample_stem(filename: str) -> str:
    """Strip all known extensions to produce a clean sample name.
    e.g. sample1.fastq.gz -> sample1
         sample2.fasta    -> sample2
    """
    for ext in SUPPORTED_EXTENSIONS:          # longest first, so .fastq.gz wins over .gz
        if filename.endswith(ext):
            return filename[: -len(ext)]
    return filename


def run_kaiju(input_file: Path, output_file: Path) -> bool:
    """Run kaiju on *input_file*, writing raw output to *output_file*.
    Returns True on success."""
    cmd = [
        KAIJU_BIN,
        "-t", NODES_FILE,
        "-f", FMI_FILE,
        "-i", str(input_file),
        "-o", str(output_file),
        "-z", str(KAIJU_THREADS),
        "-a", KAIJU_MODE,
        "-s", str(KAIJU_MIN_SCORE),
        "-e", str(KAIJU_MISMATCHES),
    ]
    print(f"  [kaiju] {input_file.name}")
    print(f"    CMD: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [ERROR] kaiju failed for {input_file.name}:\n{result.stderr}")
        return False
    return True


def run_add_taxon_names(raw_file: Path, named_file: Path) -> bool:
    """Run kaiju-addTaxonNames on *raw_file*, writing to *named_file*.
    Returns True on success."""
    cmd = [
        ADD_TAXON_NAMES_BIN,
        "-t", NODES_FILE,
        "-n", NAMES_FILE,
        "-i", str(raw_file),
        "-o", str(named_file),
    ]
    if OMIT_UNCLASSIFIED_NAMES:
        cmd.append("-u")
    if PRINT_FULL_TAXON_PATH:
        cmd.append("-p")
    elif TAXON_RANKS:
        cmd.extend(["-r", TAXON_RANKS])

    print(f"  [addTaxonNames] {raw_file.name} -> {named_file.name}")
    print(f"    CMD: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [ERROR] kaiju-addTaxonNames failed for {raw_file.name}:\n{result.stderr}")
        return False
    return True


def check_executables():
    """Abort early if required executables are not found."""
    for exe in (KAIJU_BIN, ADD_TAXON_NAMES_BIN):
        result = subprocess.run(["which", exe], capture_output=True, text=True)
        if result.returncode != 0:
            sys.exit(
                f"[ERROR] Executable not found on PATH: {exe}\n"
                f"        Set the full path in the CONFIGURATION section."
            )


def check_db_files():
    """Abort early if any database file is missing."""
    for label, path in [("FMI", FMI_FILE), ("NODES", NODES_FILE), ("NAMES", NAMES_FILE)]:
        if not Path(path).is_file():
            sys.exit(
                f"[ERROR] {label} file not found: {path}\n"
                f"        Update the path in the CONFIGURATION section."
            )


def main():
    print("=" * 60)
    print("  Kaiju Pipeline Runner")
    print("=" * 60)

    check_executables()
    check_db_files()

    # Prepare output subdirectories
    out_dir   = Path(OUTPUT_FOLDER)
    raw_dir   = out_dir / "kaiju_raw"
    named_dir = out_dir / "kaiju_named"
    for d in (out_dir, raw_dir, named_dir):
        d.mkdir(parents=True, exist_ok=True)

    # Discover input files
    input_files = find_input_files(INPUT_FOLDER)
    if not input_files:
        sys.exit(f"[ERROR] No FASTA/FASTQ/GZ files found in: {INPUT_FOLDER}")

    print(f"\nFound {len(input_files)} input file(s) in: {INPUT_FOLDER}\n")

    ok_count = fail_count = 0

    for input_file in input_files:
        stem      = sample_stem(input_file.name)
        raw_out   = raw_dir   / f"{stem}.kaiju.out"
        named_out = named_dir / f"{stem}.kaiju_names.tsv"

        print(f"\n{'─' * 50}")
        print(f"  Sample : {stem}")
        print(f"{'─' * 50}")

        if not run_kaiju(input_file, raw_out):
            fail_count += 1
            continue

        if run_add_taxon_names(raw_out, named_out):
            ok_count += 1
        else:
            fail_count += 1

    print(f"\n{'=' * 60}")
    print(f"  Done.  Successful: {ok_count}  |  Failed: {fail_count}")
    print(f"  Raw outputs  : {raw_dir}")
    print(f"  Named outputs: {named_dir}")
    print(f"{'=' * 60}\n")


if __name__ == "__main__":
    main()