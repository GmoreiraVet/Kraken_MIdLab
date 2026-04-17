#!/usr/bin/env python3

import os
import re
from Bio import SeqIO

# =========================
# CONFIGURATION SECTION
# =========================

CONFIG = {
    # Input multi-FASTA
    "input_fasta": "/home/viroicbas2023/Documents/Gmoreira/EnteroVirusCarraçasUropigeaCustom/Enterovirus/ABRIL_2026/Seqs/Database/Varios.fasta",

    # Output directory
    "output_dir": "/home/viroicbas2023/Documents/Gmoreira/EnteroVirusCarraçasUropigeaCustom/Enterovirus/ABRIL_2026/Seqs/Database/Varios_Individualizados",
    # Overwrite existing files
    "overwrite": True,

    # Remove problematic filename characters
    "sanitize_filenames": True,

    # Max filename length (avoid OS issues)
    "max_name_length": 100,
}

# =========================
# HELPERS
# =========================

def sanitize_filename(name):
    """Convert header to safe filename."""
    # Replace spaces with underscores
    name = name.replace(" ", "_")

    # Remove problematic characters
    name = re.sub(r"[^\w\-\.]", "", name)

    # Trim length
    name = name[:CONFIG["max_name_length"]]

    return name


# =========================
# MAIN LOGIC
# =========================

def main():
    input_fasta = CONFIG["input_fasta"]
    output_dir = CONFIG["output_dir"]

    os.makedirs(output_dir, exist_ok=True)

    count = 0

    for record in SeqIO.parse(input_fasta, "fasta"):

        header = record.description  # full header line

        if CONFIG["sanitize_filenames"]:
            filename = sanitize_filename(header)
        else:
            filename = header.replace(" ", "_")

        output_path = os.path.join(output_dir, f"{filename}.fasta")

        if os.path.exists(output_path) and not CONFIG["overwrite"]:
            print(f"Skipping existing file: {output_path}")
            continue

        SeqIO.write(record, output_path, "fasta")
        count += 1

    print(f"\nDone.")
    print(f"Sequences written: {count}")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()