import os
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import logging
from collections import defaultdict, Counter

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
# DIAMOND pre-filters by e-value only (--min-identity requires v2.1.0+).
# Identity, alignment length, and e-value are all enforced strictly in parse().
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
        "qstart qend sstart send evalue bitscore qseq_translated "
        "--very-sensitive "
        "--max-target-seqs 25 "
        f"--evalue {MAX_EVALUE} "   # e-value pre-filter (compatible with all versions)
        "--threads 4"
    )

    run(cmd)
    return out_tsv

# =============================================================================
# PARSE COVERAGE  (with strict post-filtering)
# TSV columns (0-based):
#   0  qseqid   1  sseqid   2  pident   3  length   4  mismatch
#   5  gapopen  6  qstart   7  qend     8  sstart   9  send
#  10  evalue  11  bitscore  12  qseq
#
# qseq is the aligned query segment already translated to AA by DIAMOND.
# Gaps in the alignment ('-') are skipped when voting.
# votes[protein][ref_pos] = Counter({'A': 5, 'V': 1, ...})
# =============================================================================

def parse(tsv, ref_lengths):
    """Build per-protein coverage arrays and AA vote tables."""

    FALLBACK = 5000

    cov   = {pid: np.zeros(length + 1, dtype=int)
             for pid, length in ref_lengths.items()}

    # votes[protein] is a list of Counters, one per reference AA position
    votes = {pid: [Counter() for _ in range(length + 1)]
             for pid, length in ref_lengths.items()}

    counters = defaultdict(int)
    unknown_proteins = set()

    total = filtered_identity = filtered_evalue = filtered_length = accepted = 0

    with open(tsv) as f:
        for line in f:

            cols = line.strip().split("\t")
            if len(cols) < 13:
                continue

            total += 1

            protein = cols[1]

            try:
                pident  = float(cols[2])
                aln_len = int(cols[3])
                sstart  = int(cols[8])   # 1-based AA position on the reference
                send    = int(cols[9])
                evalue  = float(cols[10])
            except ValueError:
                continue

            qseq = cols[12]   # AA sequence of the aligned query segment

            # ── Post-filters ─────────────────────────────────────────────────
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

            # Ensure data structures exist for unexpected proteins
            if protein not in cov:
                if protein not in unknown_proteins:
                    logging.warning(
                        f"Protein '{protein}' not in reference FASTA; "
                        f"using fallback array of {FALLBACK} aa"
                    )
                    unknown_proteins.add(protein)
                cov[protein]   = np.zeros(FALLBACK, dtype=int)
                votes[protein] = [Counter() for _ in range(FALLBACK)]

            # ── Coverage depth ────────────────────────────────────────────────
            s, e = min(sstart, send), max(sstart, send)
            arr_len = len(cov[protein])
            s = max(0, min(s, arr_len - 1))
            e = max(0, min(e, arr_len))
            cov[protein][s:e] += 1

            # ── AA vote: walk qseq char by char, skip gap characters ──────────
            ref_pos = min(sstart, send) - 1   # convert to 0-based
            for aa in qseq:
                if aa == '-':
                    ref_pos += 1              # gap in query: advance ref, no vote
                    continue
                if 0 <= ref_pos < len(votes[protein]):
                    votes[protein][ref_pos][aa] += 1
                ref_pos += 1

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

    return cov, votes


# =============================================================================
# BUILD CONSENSUS SEQUENCES
# =============================================================================

def build_consensus(votes, ref_lengths):
    """
    For each protein, build a consensus AA sequence:
      - Most common amino acid at positions with coverage
      - 'N' at positions with no coverage (gap)
      - Stop codons (*) are replaced with X to keep BLAST-friendly output
    Writes one FASTA per protein to OUTPUT_FOLDER/consensus/.
    Returns {protein: consensus_string}.
    """

    consensus_dir = os.path.join(OUTPUT_FOLDER, "consensus")
    os.makedirs(consensus_dir, exist_ok=True)

    all_consensus = {}

    for protein, pos_votes in votes.items():

        true_len = ref_lengths.get(protein, len(pos_votes))
        seq = []

        covered = 0
        for pos in range(true_len):
            counter = pos_votes[pos] if pos < len(pos_votes) else Counter()
            if not counter:
                seq.append('N')
            else:
                best = counter.most_common(1)[0][0]
                best = 'X' if best == '*' else best   # no stop codons in BLAST query
                seq.append(best)
                covered += 1

        consensus = "".join(seq)
        all_consensus[protein] = consensus

        pct = 100.0 * covered / true_len if true_len else 0
        logging.info(
            f"Consensus {protein}: {covered}/{true_len} positions covered "
            f"({pct:.1f}%), {consensus.count('N')} N gaps"
        )

        # Write FASTA — wrap at 60 chars per line
        safe_name = protein.replace("/", "_").replace(" ", "_")
        out_path = os.path.join(consensus_dir, f"{safe_name}_consensus.faa")

        with open(out_path, "w") as fh:
            fh.write(f">{protein}_consensus  gaps=N  filters=id{MIN_IDENTITY}_e{MAX_EVALUE}_aln{MIN_ALN_LEN}\n")
            for i in range(0, len(consensus), 60):
                fh.write(consensus[i:i+60] + "\n")

        logging.info(f"Consensus FASTA saved: {out_path}")

    return all_consensus

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
 
    cov, votes = parse(tsv, ref_lengths)
 
    plot(cov, ref_lengths)
 
    build_consensus(votes, ref_lengths)
 
    logging.info("Pipeline complete")
 
    print("""
　　　　　   __
　　　　 ／フ   フ
　　　　|  .   .|
　 　　／`ミ__xノ
　 　 /　　 　 |
　　 /　 ヽ　　ﾉ
 　 │　　 | | |
／￣|　　 | | |
| (￣ヽ_ヽ)_)__)
＼二つ
 
  Mimi did her best to assemble your genome.
""")
 
# =============================================================================
# ENTRY
# =============================================================================
 
if __name__ == "__main__":
    main()