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
# QUALITY FILTERS
# Adjust these to control stringency:
#   MIN_IDENTITY  : minimum % amino acid identity (0–100). 70 = same viral species.
#   MAX_EVALUE    : maximum e-value. 1e-5 is tighter than DIAMOND's default 1e-3.
#   MIN_ALN_LEN   : minimum alignment length in AA. Removes tiny spurious hits.
# =============================================================================

MIN_IDENTITY = 70.0   # %
MAX_EVALUE   = 1e-5
MIN_ALN_LEN  = 25     # amino acids

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
    run(f"diamond makedb --in {REFERENCE_FASTA} -d {DB_PREFIX}")

    if not os.path.exists(dmnd_file):
        raise RuntimeError("DIAMOND DB creation failed")

    logging.info("Database built successfully")
    return dmnd_file

# =============================================================================
# READ REFERENCE LENGTHS FROM FASTA
# =============================================================================

def get_reference_lengths(fasta_path):
    """Return {seq_id: length_in_aa} for every entry in the reference FASTA."""

    lengths = {}
    current_id = None
    current_len = 0

    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_len
                current_id = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line)

    if current_id is not None:
        lengths[current_id] = current_len

    logging.info(f"Reference proteins loaded: {len(lengths)}")
    for pid, l in lengths.items():
        logging.info(f"  {pid}: {l} aa")

    return lengths

# =============================================================================
# FASTQ → FASTA
# =============================================================================

def fastq_to_fasta(fastq, fasta_out):

    run(f"seqtk seq -a {fastq} > {fasta_out}")
    return fasta_out

# =============================================================================
# RUN DIAMOND
# Post-filtering in parse() is the strict second pass.
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
        f"--evalue {MAX_EVALUE} "          # tighter e-value pre-filter
        "--threads 4"
    )

    run(cmd)
    return out_tsv

# =============================================================================
# PARSE COVERAGE  (with strict post-filtering)
# TSV columns (0-based):
#   0  qseqid   1  sseqid   2  pident   3  length   4  mismatch
#   5  gapopen  6  qstart   7  qend     8  sstart   9  send
#  10  evalue  11  bitscore
# =============================================================================

def parse(tsv, ref_lengths):
    """Build per-protein coverage arrays, keeping only high-confidence hits."""

    FALLBACK = 5000

    cov = {pid: np.zeros(length + 1, dtype=int)
           for pid, length in ref_lengths.items()}

    counters = defaultdict(int)   # per-protein accepted hit counter
    unknown_proteins = set()

    total = filtered_identity = filtered_evalue = filtered_length = accepted = 0

    with open(tsv) as f:
        for line in f:

            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue

            total += 1

            protein = cols[1]

            try:
                pident  = float(cols[2])
                aln_len = int(cols[3])
                start   = int(cols[8])
                end     = int(cols[9])
                evalue  = float(cols[10])
            except ValueError:
                continue

            # ── Post-filters (belt-and-braces after DIAMOND pre-filters) ─────
            if pident < MIN_IDENTITY:
                filtered_identity += 1
                continue

            if evalue > MAX_EVALUE:
                filtered_evalue += 1
                continue

            if aln_len < MIN_ALN_LEN:
                filtered_length += 1
                continue
            # ─────────────────────────────────────────────────────────────────

            accepted += 1
            counters[protein] += 1

            if protein not in cov:
                if protein not in unknown_proteins:
                    logging.warning(
                        f"Protein '{protein}' not in reference FASTA; "
                        f"using fallback array of {FALLBACK} aa"
                    )
                    unknown_proteins.add(protein)
                cov[protein] = np.zeros(FALLBACK, dtype=int)

            s, e = min(start, end), max(start, end)
            arr_len = len(cov[protein])
            s = max(0, min(s, arr_len - 1))
            e = max(0, min(e, arr_len))
            cov[protein][s:e] += 1

    # ── Summary log ──────────────────────────────────────────────────────────
    logging.info("─── Hit filtering summary ───────────────────────────────")
    logging.info(f"  Total alignments parsed  : {total}")
    logging.info(f"  Removed (pident < {MIN_IDENTITY}%)  : {filtered_identity}")
    logging.info(f"  Removed (evalue > {MAX_EVALUE})  : {filtered_evalue}")
    logging.info(f"  Removed (aln_len < {MIN_ALN_LEN} aa): {filtered_length}")
    logging.info(f"  Accepted (high-confidence): {accepted}")
    logging.info("  Per-protein accepted hits:")
    for pid, n in sorted(counters.items()):
        logging.info(f"    {pid}: {n} hits")
    logging.info("─────────────────────────────────────────────────────────")

    return cov

# =============================================================================
# SMOOTH
# =============================================================================

def smooth(x, w=20):

    if len(x) < w:
        return x

    return np.convolve(x, np.ones(w) / w, mode="valid")

# =============================================================================
# PLOT
# =============================================================================

def plot(cov_dict, ref_lengths):

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    fill_color = line_color = CUSTOM_COLOR or "#2A9D8F"

    if not any(c is not None and np.sum(c) > 0 for c in cov_dict.values()):

        plt.figure(figsize=(10, 5))
        plt.text(0.5, 0.5, "No high-confidence DIAMOND hits",
                 ha="center", va="center", fontsize=16)
        out = os.path.join(OUTPUT_FOLDER, "no_hits.png")
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()
        logging.warning("No high-confidence hits to plot")
        return

    for protein, c in cov_dict.items():

        if c is None or np.sum(c) == 0:
            logging.info(f"Skipping {protein}: zero coverage after filtering")
            continue

        true_len  = ref_lengths.get(protein, len(c))
        c_trimmed = c[:true_len]

        sm = smooth(c_trimmed)
        x  = np.arange(len(sm))

        fig, ax = plt.subplots(figsize=(12, 4))

        ax.plot(x, sm, color=line_color, linewidth=1.2)
        ax.fill_between(x, sm, color=fill_color, alpha=0.4)

        ax.set_xlim(0, true_len)
        ax.set_ylim(bottom=0)

        ax.set_title(f"Protein coverage: {protein}")
        ax.set_xlabel("Amino acid position")
        ax.set_ylabel("Read depth")

        # Filters used — shown on the plot so figures are self-documenting
        filter_label = (
            f"Filters: identity ≥ {MIN_IDENTITY}%  |  "
            f"e-value ≤ {MAX_EVALUE}  |  "
            f"aln ≥ {MIN_ALN_LEN} aa  |  "
            f"ref length: {true_len} aa"
        )
        ax.text(
            0.5, 1.02, filter_label,
            transform=ax.transAxes,
            ha="center", va="bottom",
            fontsize=7, color="grey"
        )

        safe_name = protein.replace("/", "_").replace(" ", "_")
        out = os.path.join(OUTPUT_FOLDER, f"{safe_name}_pileup.png")
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()

        logging.info(f"Saved {out}")

# =============================================================================
# MAIN
# =============================================================================

def main():

    setup_logging()

    logging.info(
        f"Quality filters — identity: ≥{MIN_IDENTITY}%  "
        f"e-value: ≤{MAX_EVALUE}  aln_len: ≥{MIN_ALN_LEN} aa"
    )

    ref_lengths = get_reference_lengths(REFERENCE_FASTA)

    dmnd = build_db_if_needed()

    tsv = run_diamond(dmnd)

    cov = parse(tsv, ref_lengths)

    plot(cov, ref_lengths)

    logging.info("Pipeline complete")

# =============================================================================
# ENTRY
# =============================================================================

if __name__ == "__main__":
    main()