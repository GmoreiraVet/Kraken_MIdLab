#!/usr/bin/env python3

import os
import shutil
import subprocess
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

# Input fasta/contigs file
INPUT_FASTA = "/home/viroicbas2023/Documents/Gmoreira/Carraças_Metagenomica_Gui/Reads_DeRNAfied/VirusNOVOS/42/PhlebovirusNovo/PVGA/Resultados/L_Segment/pvga_42_filtered.fastq_ON_TickPhleboVirus_L_Segment.fa"

# Protein reference database
REFERENCE_PROTEINS = "/home/viroicbas2023/Documents/Gmoreira/Carraças_Metagenomica_Gui/Reads_DeRNAfied/VirusNOVOS/42/PhlebovirusNovo/ProovFrame/BrownDogTickPhlebovirus-Reference.faa"

# Working directory for intermediate files
WORK_DIR = "/home/viroicbas2023/Documents/Gmoreira/Carraças_Metagenomica_Gui/Reads_DeRNAfied/VirusNOVOS/42/PhlebovirusNovo/ProovFrame"

# Final output directory
FINAL_OUTPUT_DIR = "/home/viroicbas2023/Documents/Gmoreira/CorrectedGenomes"

# Number of threads
THREADS = 30

# Intermediate mapping output
MAP_OUTPUT = "PhleboRef.tsv"

# Corrected fasta output
CORRECTED_FASTA = "CorrectedPhlebo_LSegment.fasta"

# ============================================================
# FUNCTIONS
# ============================================================

def run_command(cmd):
    """Run a shell command and stop on failure."""
    print("\nRunning:")
    print(" ".join(cmd))
    print()

    result = subprocess.run(cmd)

    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {result.returncode}"
        )

# ============================================================
# MAIN
# ============================================================

def main():

    os.makedirs(WORK_DIR, exist_ok=True)
    os.makedirs(FINAL_OUTPUT_DIR, exist_ok=True)

    map_output_path = os.path.join(WORK_DIR, MAP_OUTPUT)
    corrected_output_path = os.path.join(WORK_DIR, CORRECTED_FASTA)

    # --------------------------------------------------------
    # Step 1: proovframe map
    # --------------------------------------------------------

    map_cmd = [
        "proovframe",
        "map",
        "-a", REFERENCE_PROTEINS,
        "-o", map_output_path,
        "-t", str(THREADS),
        INPUT_FASTA
    ]

    run_command(map_cmd)

    # --------------------------------------------------------
    # Step 2: proovframe fix
    # --------------------------------------------------------

    fix_cmd = [
        "proovframe",
        "fix",
        "-o", corrected_output_path,
        INPUT_FASTA,
        map_output_path
    ]

    run_command(fix_cmd)

    # --------------------------------------------------------
    # Copy final result
    # --------------------------------------------------------

    final_destination = os.path.join(
        FINAL_OUTPUT_DIR,
        os.path.basename(corrected_output_path)
    )

    shutil.copy2(corrected_output_path, final_destination)

    print("\nFinished successfully.")
    print(f"Final output:")
    print(final_destination)


if __name__ == "__main__":
    main()