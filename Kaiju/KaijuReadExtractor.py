#!/usr/bin/env python3

import gzip

# =========================================================
# 🛠️ USER SETTINGS
# =========================================================

KAIJU_FILE = "/home/viroicbas2023/Documents/Gmoreira/AndreiaGuiMay8/KaijuResults_BIGDB/kaiju_named/KBFTRN_12_sample_12.kaiju_names.tsv"
FASTQ_FILE = "/home/viroicbas2023/Documents/Gmoreira/AndreiaGuiMay8/KBFTRN_fastq/KBFTRN_12_sample_12.fastq.gz"

# Output file
OUTPUT_FILE = "/home/viroicbas2023/Downloads/Fotos/HPV_GUI.fasta"

# Add ONE or MULTIPLE taxIDs
TARGET_TAXIDS = {
    "10566",
}

# Output mode:
# "fasta" or "fastq"
OUTPUT_FORMAT = "fasta"

KAIJU_IS_GZIPPED = False
FASTQ_IS_GZIPPED = True

# =========================================================


def open_file(path, gz=False):
    if gz or path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def extract_read_ids():
    """
    Extract read IDs from Kaiju extended output format:
    C <readID> <taxID> ...
    """
    read_ids = set()

    with open_file(KAIJU_FILE, KAIJU_IS_GZIPPED) as f:
        for line in f:
            if not line.strip():
                continue

            cols = line.rstrip().split("\t")

            # Skip unclassified reads
            if cols[0] != "C":
                continue

            if len(cols) < 3:
                continue

            read_id = cols[1]
            taxid = cols[2]

            # Match any target taxID
            if taxid in TARGET_TAXIDS:
                read_ids.add(read_id)

    return read_ids


def extract_reads(read_ids):
    """
    Extract matching reads from FASTQ.
    Output can be FASTA or FASTQ.
    """

    with open_file(FASTQ_FILE, FASTQ_IS_GZIPPED) as fq, open(OUTPUT_FILE, "w") as out:

        while True:
            header = fq.readline().strip()

            if not header:
                break

            seq = fq.readline().strip()
            plus = fq.readline().strip()
            qual = fq.readline().strip()

            read_id = header.split()[0][1:]  # remove '@'

            if read_id in read_ids:

                if OUTPUT_FORMAT.lower() == "fasta":
                    out.write(f">{read_id}\n{seq}\n")

                elif OUTPUT_FORMAT.lower() == "fastq":
                    out.write(f"{header}\n{seq}\n{plus}\n{qual}\n")

                else:
                    raise ValueError(
                        "OUTPUT_FORMAT must be 'fasta' or 'fastq'"
                    )


def main():

    print("🔍 Parsing Kaiju output...")
    read_ids = extract_read_ids()

    print(
        f"✔ Found {len(read_ids)} reads for "
        f"{len(TARGET_TAXIDS)} target taxIDs"
    )

    print("✂ Extracting reads...")
    extract_reads(read_ids)

    print(f"✅ Done → {OUTPUT_FILE}")


if __name__ == "__main__":
    main()