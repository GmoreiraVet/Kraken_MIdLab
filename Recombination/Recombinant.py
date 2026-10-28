#!/usr/bin/env python3
"""
recomb_split_trees.py
=====================
For each user-defined alignment interval (putative recombinant region):

  1. Extract the interval          -> recombinant.fasta
  2. Extract everything outside it -> non_recombinant.fasta
     (left flank + right flank, concatenated)
  3. Run IQ-TREE 3 on both         -> iqtree3 -s <file>.phy --alrt 10000 -B 10000

Output layout
-------------
OUTPUT_DIR/
├── summary.tsv
├── run.log
└── interval_01_<start>-<end>[_<label>]/
    ├── recombinant/
    │   ├── recombinant.fasta
    │   ├── recombinant.phy
    │   └── recombinant.phy.treefile (+ other IQ-TREE files)
    └── non_recombinant/
        ├── non_recombinant.fasta
        ├── non_recombinant.phy
        └── non_recombinant.phy.treefile (+ other IQ-TREE files)

Coordinates are 1-based, inclusive, and refer to ALIGNMENT columns
(i.e. including gaps), as displayed in e.g. AliView / Geneious / MEGA.

Requirements: Python >= 3.8, iqtree3 in PATH (no third-party Python packages).
Usage: edit the CONFIG block below, then run the script (F5 in VS Code).
"""

# =============================================================================
# CONFIG  — edit everything here
# =============================================================================
from pathlib import Path

# Input alignment (aligned FASTA, all sequences same length)
ALIGNMENT_FILE = Path("MolecularClock_Analysis/L_Segment/Recombination_Analysis/L_Segment_Database_ALIGNED_WithPTSequence_TRIIMED_FOR_RDP.fasta")

# Parent output folder (created if missing)
OUTPUT_DIR = Path("/home/viroicbas2023/Documents/Gmoreira/Carraças/Carraças_Metagenomica_Gui/VOLZSHKOESEPARADO/MolecularClock_Analysis/L_Segment/Recombination_Analysis/Trees")

# Intervals to analyse: (start, end) or (start, end, "label"). 1-based, inclusive.
INTERVALS = [
    (148, 1034, "PP473448"),
    (6071, 6930, "PP069765"),
    (2663, 4120, "PP069785")
    # (8000, 9100, "region_C"),
]

# --- IQ-TREE ---------------------------------------------------------------
RUN_IQTREE = True                  # False = only write FASTA/PHYLIP files
IQTREE_EXECUTABLE = "iqtree3"      # name in PATH, or full path
IQTREE_ARGS = ["--alrt", "10000", "-B", "10000"]
IQTREE_THREADS = "AUTO"            # "AUTO" or an integer (as string), passed to -T
IQTREE_EXTRA_ARGS = []             # e.g. ["-m", "GTR+G", "-o", "OutgroupName"]
TREE_INPUT_FORMAT = "fasta"       # "phylip" (writes .phy, as in your example) or "fasta"
SKIP_IF_TREE_EXISTS = True         # skip IQ-TREE if a .treefile is already present

# --- Sequence handling -----------------------------------------------------
# Drop sequences with no informative sites (only gaps/N/?) in a given region.
DROP_EMPTY_SEQUENCES = True
# Warn if a region is shorter than this (alignment columns)
MIN_REGION_LENGTH = 100
# FASTA line wrapping (0 = single-line sequences)
FASTA_LINE_WIDTH = 80
# =============================================================================


import shutil
import subprocess
import sys
from datetime import datetime

LOG_FILE = None  # set in main()


def log(msg: str) -> None:
    line = f"[{datetime.now():%H:%M:%S}] {msg}"
    print(line, flush=True)
    if LOG_FILE is not None:
        with open(LOG_FILE, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")


# ----------------------------------------------------------------------------
# I/O helpers
# ----------------------------------------------------------------------------
def read_fasta(path: Path) -> dict:
    """Read FASTA into an ordered dict {id: sequence}. ID = first word of header."""
    seqs, name, buf = {}, None, []
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0] if line[1:].strip() else ""
                if not name:
                    sys.exit("ERROR: empty FASTA header found.")
                if name in seqs:
                    sys.exit(f"ERROR: duplicate sequence ID '{name}'.")
                buf = []
            else:
                if name is None:
                    sys.exit("ERROR: sequence data before first FASTA header.")
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def write_fasta(seqs: dict, path: Path, width: int = 80) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            if width and width > 0:
                for i in range(0, len(seq), width):
                    fh.write(seq[i:i + width] + "\n")
            else:
                fh.write(seq + "\n")


def write_phylip(seqs: dict, path: Path) -> None:
    """Relaxed sequential PHYLIP (names not truncated; IQ-TREE reads this)."""
    n = len(seqs)
    length = len(next(iter(seqs.values())))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f" {n} {length}\n")
        for name, seq in seqs.items():
            fh.write(f"{name} {seq}\n")


def is_empty(seq: str) -> bool:
    return all(c in "-?Nn" for c in seq)


# ----------------------------------------------------------------------------
# Core logic
# ----------------------------------------------------------------------------
def validate_intervals(intervals, aln_len):
    parsed = []
    for i, iv in enumerate(intervals, 1):
        if len(iv) == 2:
            start, end = iv
            label = ""
        elif len(iv) == 3:
            start, end, label = iv
        else:
            sys.exit(f"ERROR: interval #{i} must be (start, end) or (start, end, label).")
        if not (isinstance(start, int) and isinstance(end, int)):
            sys.exit(f"ERROR: interval #{i}: start/end must be integers.")
        if start < 1 or end > aln_len or start > end:
            sys.exit(f"ERROR: interval #{i} ({start}-{end}) invalid for alignment length {aln_len}.")
        parsed.append((start, end, str(label).strip().replace(" ", "_")))
    return parsed


def split_alignment(seqs: dict, start: int, end: int):
    """Return (recombinant, non_recombinant) dicts. start/end 1-based inclusive."""
    s0, e0 = start - 1, end  # python slice
    rec = {k: v[s0:e0] for k, v in seqs.items()}
    non = {k: v[:s0] + v[e0:] for k, v in seqs.items()}
    return rec, non


def prepare_region(seqs: dict, tag: str, outdir: Path):
    """Write FASTA (+PHYLIP) for one region; return path to file for IQ-TREE or None."""
    outdir.mkdir(parents=True, exist_ok=True)
    kept = seqs
    if DROP_EMPTY_SEQUENCES:
        kept = {k: v for k, v in seqs.items() if not is_empty(v)}
        dropped = [k for k in seqs if k not in kept]
        if dropped:
            log(f"  [{tag}] dropped {len(dropped)} empty sequence(s): {', '.join(dropped)}")

    length = len(next(iter(kept.values()))) if kept else 0
    if len(kept) < 4:
        log(f"  [{tag}] WARNING: only {len(kept)} sequences — IQ-TREE needs >= 4. Skipping tree.")
        if kept:
            write_fasta(kept, outdir / f"{tag}.fasta", FASTA_LINE_WIDTH)
        return None, len(kept), length
    if length < MIN_REGION_LENGTH:
        log(f"  [{tag}] WARNING: region is only {length} columns (< {MIN_REGION_LENGTH}).")

    fasta = outdir / f"{tag}.fasta"
    write_fasta(kept, fasta, FASTA_LINE_WIDTH)
    if TREE_INPUT_FORMAT.lower() == "phylip":
        phy = outdir / f"{tag}.phy"
        write_phylip(kept, phy)
        return phy, len(kept), length
    return fasta, len(kept), length


def run_iqtree(infile: Path) -> bool:
    cmd = [IQTREE_EXECUTABLE, "-s", infile.name, *IQTREE_ARGS,
           "-T", str(IQTREE_THREADS), *IQTREE_EXTRA_ARGS]
    treefile = infile.parent / (infile.name + ".treefile")
    if SKIP_IF_TREE_EXISTS and treefile.exists():
        log(f"  treefile exists, skipping: {treefile}")
        return True
    log(f"  running: {' '.join(cmd)}  (cwd={infile.parent})")
    try:
        # cwd = region folder so all IQ-TREE outputs stay there
        res = subprocess.run(cmd, cwd=infile.parent, capture_output=True, text=True)
    except FileNotFoundError:
        log(f"  ERROR: '{IQTREE_EXECUTABLE}' not found.")
        return False
    (infile.parent / "iqtree_stdout.log").write_text(res.stdout + "\n" + res.stderr, encoding="utf-8")
    if res.returncode != 0:
        log(f"  ERROR: IQ-TREE exited with code {res.returncode} "
            f"(see {infile.parent / 'iqtree_stdout.log'})")
        return False
    log(f"  done: {treefile}")
    return True


def main() -> None:
    global LOG_FILE
    if not ALIGNMENT_FILE.is_file():
        sys.exit(f"ERROR: alignment not found: {ALIGNMENT_FILE.resolve()}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_FILE = OUTPUT_DIR / "run.log"

    if RUN_IQTREE and shutil.which(IQTREE_EXECUTABLE) is None:
        sys.exit(f"ERROR: '{IQTREE_EXECUTABLE}' not found in PATH. "
                 f"Set IQTREE_EXECUTABLE or RUN_IQTREE = False.")

    seqs = read_fasta(ALIGNMENT_FILE)
    if not seqs:
        sys.exit("ERROR: no sequences read.")
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        sys.exit(f"ERROR: sequences have different lengths ({sorted(lengths)[:5]}...) — not an alignment.")
    aln_len = lengths.pop()
    log(f"Alignment: {len(seqs)} sequences x {aln_len} columns")

    intervals = validate_intervals(INTERVALS, aln_len)
    summary = ["interval\tstart\tend\tlabel\trec_len\tnonrec_len\trec_nseq\tnonrec_nseq\trec_tree\tnonrec_tree"]

    for idx, (start, end, label) in enumerate(intervals, 1):
        name = f"interval_{idx:02d}_{start}-{end}" + (f"_{label}" if label else "")
        base = OUTPUT_DIR / name
        log(f"=== {name} ===")
        rec, non = split_alignment(seqs, start, end)

        results = {}
        for tag, data in (("recombinant", rec), ("non_recombinant", non)):
            infile, n, length = prepare_region(data, tag, base / tag)
            ok = None
            if infile is not None and RUN_IQTREE:
                ok = run_iqtree(infile)
            results[tag] = (n, length, "NA" if ok is None else ("ok" if ok else "FAILED"))

        summary.append("\t".join(map(str, [
            name, start, end, label or "-",
            results["recombinant"][1], results["non_recombinant"][1],
            results["recombinant"][0], results["non_recombinant"][0],
            results["recombinant"][2], results["non_recombinant"][2],
        ])))

    (OUTPUT_DIR / "summary.tsv").write_text("\n".join(summary) + "\n", encoding="utf-8")
    log(f"Finished. Summary: {OUTPUT_DIR / 'summary.tsv'}")


if __name__ == "__main__":
    main()