import os
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import logging
from collections import defaultdict

# =============================================================================
# CONFIG
# =============================================================================

OUTPUT_FOLDER = "/home/viroicbas2023/Documents/Gmoreira/AndreiaGuiMay8/Diamond_Mimiplot_Teste"

REFERENCE_FASTA = "/home/viroicbas2023/Documents/Gmoreira/AndreiaGuiMay8/Diamond_Mimiplot_Teste/Refs/Proteins_ClosestMatch.fasta"

QUERY_FASTQ = "/home/viroicbas2023/Documents/Gmoreira/AndreiaGuiMay8/KBFTRN_fastq/KBFTRN_2_sample_02.fastq.gz"

DB_PREFIX = os.path.join(OUTPUT_FOLDER, "protein_db")

CUSTOM_COLOR = "#00FFAE"

# =============================================================================
# LOGGING
# =============================================================================

def setup_logging():

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    logging.basicConfig(
        filename=os.path.join(OUTPUT_FOLDER, "pipeline.log"),
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

    logging.getLogger().addHandler(logging.StreamHandler())

# =============================================================================
# RUN COMMAND
# =============================================================================

def run(cmd):

    logging.info(cmd)

    subprocess.run(
        cmd,
        shell=True,
        check=True,
        executable="/bin/bash"
    )

# =============================================================================
# BUILD DIAMOND DB IF NEEDED
# =============================================================================

def build_db_if_needed():

    dmnd_file = DB_PREFIX + ".dmnd"

    if os.path.exists(dmnd_file):

        logging.info("DIAMOND DB already exists, skipping build")

        return dmnd_file

    logging.info("Building DIAMOND database...")

    cmd = f"diamond makedb --in {REFERENCE_FASTA} -d {DB_PREFIX}"

    run(cmd)

    if not os.path.exists(dmnd_file):

        raise RuntimeError("DIAMOND DB creation failed")

    logging.info("Database built successfully")

    return dmnd_file

# =============================================================================
# FASTQ → FASTA
# =============================================================================

def fastq_to_fasta(fastq, fasta_out):

    run(f"seqtk seq -a {fastq} > {fasta_out}")

    return fasta_out

# =============================================================================
# RUN DIAMOND
# =============================================================================

def run_diamond(dmnd):

    fasta_query = os.path.join(OUTPUT_FOLDER, "query.fasta")

    fastq_to_fasta(QUERY_FASTQ, fasta_query)

    out_tsv = os.path.join(OUTPUT_FOLDER, "diamond.tsv")

    cmd = (
        "diamond blastx "
        f"-q {fasta_query} "
        f"-d {dmnd} "
        f"-o {out_tsv} "
        "--outfmt 6 qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore "
        "--very-sensitive "
        "--max-target-seqs 25 "
        "--evalue 1e-3 "
        "--threads 4"
    )

    run(cmd)

    return out_tsv

# =============================================================================
# PARSE COVERAGE
# =============================================================================

def parse(tsv):

    cov = defaultdict(lambda: None)

    hits = 0
    total = 0

    with open(tsv) as f:

        for line in f:

            cols = line.strip().split("\t")

            if len(cols) < 12:
                continue

            total += 1

            protein = cols[1]

            try:
                start = int(cols[8])
                end = int(cols[9])
            except:
                continue

            hits += 1

            if cov[protein] is None:
                cov[protein] = np.zeros(5000, dtype=int)

            s, e = min(start, end), max(start, end)

            cov[protein][s:e] += 1

    logging.info(f"Total alignments: {total}")
    logging.info(f"Valid hits: {hits}")

    return cov

# =============================================================================
# SMOOTH
# =============================================================================

def smooth(x, w=20):

    if len(x) < w:
        return x

    return np.convolve(x, np.ones(w)/w, mode="valid")

# =============================================================================
# PLOT
# =============================================================================

def plot(cov_dict):

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    colors = ["#2A9D8F", "#E76F51", "#264653", "#F4A261"]

    if len(cov_dict) == 0:

        plt.figure(figsize=(10,5))
        plt.text(0.5, 0.5, "No DIAMOND hits",
                 ha="center", va="center", fontsize=16)

        out = os.path.join(OUTPUT_FOLDER, "no_hits.png")
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()

        logging.warning("No hits detected")
        return

    for protein, c in cov_dict.items():

        if c is None or np.sum(c) == 0:
            continue

        sm = smooth(c)

        plt.figure(figsize=(12,4))

        plt.plot(sm, color=CUSTOM_COLOR or colors[0])
        plt.fill_between(range(len(sm)), sm, alpha=0.4)

        plt.title(f"Protein coverage: {protein}")
        plt.xlabel("Amino acid position")
        plt.ylabel("Read depth")

        out = os.path.join(
            OUTPUT_FOLDER,
            f"{protein}_pileup.png"
        )

        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()

        logging.info(f"Saved {out}")

# =============================================================================
# MAIN
# =============================================================================

def main():

    setup_logging()

    dmnd = build_db_if_needed()

    tsv = run_diamond(dmnd)

    cov = parse(tsv)

    plot(cov)

    logging.info("Pipeline complete")

# =============================================================================
# ENTRY
# =============================================================================

if __name__ == "__main__":
    main()