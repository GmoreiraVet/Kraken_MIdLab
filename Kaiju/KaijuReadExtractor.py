#!/usr/bin/env python3

import gzip

# =========================================================
# 🛠️  USER SETTINGS
# =========================================================

KAIJU_FILE  = "/home/viroicbas2023/Documents/Gmoreira/Carraças_Metagenomica_Gui/Reads_DeRNAfied/Flye_Denovo_No_rRna/Assemblies_collected/Kaiju/kaiju_named/assembly_86.kaiju_names.tsv"
INPUT_FILE  = "/home/viroicbas2023/Documents/Gmoreira/Carraças_Metagenomica_Gui/Reads_DeRNAfied/Flye_Denovo_No_rRna/Assemblies_collected/assembly_86.fasta"
OUTPUT_FILE = "/home/viroicbas2023/Documents/Gmoreira/Carraças_Metagenomica_Gui/Reads_DeRNAfied/VirusNOVOS/86/Naiorvirus/ExtraidosKaiju_Flye_86_Naiorvirus.fasta"

TARGET_TAXIDS = {
    "2304647",
    "1971396",
    "1923252",
    "1923253",
    "1923254",
    "1980415",
    "1526514",
    "1608083",
    "1608084",
    "1608090",
}

# Output mode: "fasta" or "fastq"
OUTPUT_FORMAT    = "fasta"
KAIJU_IS_GZIPPED = False
INPUT_IS_GZIPPED = False

# =========================================================


def open_file(path, gz=False):
    if gz or path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def detect_format(path):
    """Auto-detect FASTA or FASTQ by peeking at the first non-empty line."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                return "fasta"
            if line.startswith("@"):
                return "fastq"
            raise ValueError(f"Cannot detect format — unexpected first character: {line[:30]}")
    raise ValueError("Input file appears to be empty.")


def extract_read_ids():
    """Parse Kaiju TSV and return the set of read IDs matching TARGET_TAXIDS."""
    read_ids = set()
    with open_file(KAIJU_FILE, KAIJU_IS_GZIPPED) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if cols[0] != "C" or len(cols) < 3:
                continue
            if cols[2] in TARGET_TAXIDS:
                read_ids.add(cols[1])
    return read_ids


def extract_fasta(read_ids):
    """Extract reads from FASTA — handles multi-line sequences."""
    matched, total = 0, 0
    with open_file(INPUT_FILE, INPUT_IS_GZIPPED) as fh, open(OUTPUT_FILE, "w") as out:
        current_id   = None
        seq_lines    = []
        write_record = False

        for raw_line in fh:
            line = raw_line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    total += 1
                    if write_record:
                        out.write(f">{current_id}\n{''.join(seq_lines)}\n")
                        matched += 1
                current_id   = line[1:].split()[0]
                seq_lines    = []
                write_record = current_id in read_ids
            elif line:
                seq_lines.append(line)

        # flush last record
        if current_id is not None:
            total += 1
            if write_record:
                out.write(f">{current_id}\n{''.join(seq_lines)}\n")
                matched += 1

    return matched, total


def extract_fastq(read_ids):
    """Extract reads from FASTQ — strict 4-line records."""
    matched, total = 0, 0
    with open_file(INPUT_FILE, INPUT_IS_GZIPPED) as fh, open(OUTPUT_FILE, "w") as out:
        while True:
            header = fh.readline()
            if not header:
                break
            seq  = fh.readline().rstrip()
            plus = fh.readline().rstrip()
            qual = fh.readline().rstrip()
            header = header.rstrip()
            if not header:
                continue
            total  += 1
            read_id = header.split()[0][1:]   # strip leading '@'
            if read_id in read_ids:
                matched += 1
                if OUTPUT_FORMAT.lower() == "fasta":
                    out.write(f">{read_id}\n{seq}\n")
                else:
                    out.write(f"{header}\n{seq}\n{plus}\n{qual}\n")
    return matched, total


def main():
    print("🔍 Parsing Kaiju output...")
    read_ids = extract_read_ids()
    print(f"✔  Found {len(read_ids)} read IDs across {len(TARGET_TAXIDS)} target taxIDs")

    fmt = detect_format(INPUT_FILE)
    print(f"🧬 Input detected as: {fmt.upper()}")

    print("✂  Extracting reads...")
    if fmt == "fasta":
        matched, total = extract_fasta(read_ids)
    else:
        matched, total = extract_fastq(read_ids)

    print(f"✔  Scanned {total} records — matched {matched}")

    if matched < len(read_ids):
        print(f"⚠️  {len(read_ids) - matched} Kaiju read ID(s) not found in input.")
        print("    → Check for ID suffix differences (e.g. /1, /2).")

    print(f"✅ Done → {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
