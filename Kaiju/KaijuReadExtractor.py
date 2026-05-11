#!/usr/bin/env python3

import gzip

# =========================================================
# 🛠️ USER SETTINGS
# =========================================================

KAIJU_FILE = "kaiju.out"
FASTQ_FILE = "reads.fastq"
OUTPUT_FASTA = "output.fasta"

TARGET_TAXID = "660955"   # change this to your target

KAIJU_IS_GZIPPED = False
FASTQ_IS_GZIPPED = False

# =========================================================


def open_file(path, gz=False):
    if gz or path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def extract_read_ids():
    """
    Extract read IDs from Kaiju extended format:
    C <readID> <taxID> ...
    """
    read_ids = set()

    with open_file(KAIJU_FILE, KAIJU_IS_GZIPPED) as f:
        for line in f:
            if not line.strip():
                continue

            cols = line.rstrip().split("\t")

            # Skip non-classified lines if any
            if cols[0] != "C":
                continue

            if len(cols) < 3:
                continue

            read_id = cols[1]
            taxid = cols[2]

            if taxid == TARGET_TAXID:
                read_ids.add(read_id)

    return read_ids


def fastq_to_fasta(read_ids):
    """
    Convert matching FASTQ reads to FASTA.
    """
    with open_file(FASTQ_FILE, FASTQ_IS_GZIPPED) as fq, open(OUTPUT_FASTA, "w") as out:

        while True:
            header = fq.readline().strip()
            if not header:
                break

            seq = fq.readline().strip()
            fq.readline()  # +
            fq.readline()  # quality

            read_id = header.split()[0][1:]  # remove '@'

            if read_id in read_ids:
                out.write(f">{read_id}\n{seq}\n")


def main():
    print("🔍 Parsing Kaiju extended output...")
    read_ids = extract_read_ids()

    print(f"✔ Found {len(read_ids)} reads for taxID {TARGET_TAXID}")

    print("✂ Extracting FASTQ reads...")
    fastq_to_fasta(read_ids)

    print(f"✅ Done → {OUTPUT_FASTA}")


if __name__ == "__main__":
    main()