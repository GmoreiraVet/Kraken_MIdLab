#!/usr/bin/env python3
"""
Fasta_Gene_Excisor_v3.py

Pipeline:
1. Parse GenBank → gene coordinates
2. Extract gene regions from alignment (reference projected via MAFFT --add)
3. MAFFT realignment per gene
4. IQ-TREE 3 phylogeny per gene

April 2026 — Guilherme Moreira
"""

# ==============================================================================
# CONFIGURATION
# ==============================================================================

GBK_FILE   = "/home/viroicbas2023/Downloads/Submission3073203.gbk"
FASTA_FILE = "/home/viroicbas2023/Documents/Gmoreira/EnteroVirusCarraçasUropigeaCustom/Enterovirus/ABRIL_2026/ARTIGO_NATURE_SP/Tree/Built_fixed_Aligned_Sanitized.fasta"
OUTPUT_DIR = "/home/viroicbas2023/Documents/Gmoreira/EnteroVirusCarraçasUropigeaCustom/Enterovirus/ABRIL_2026/ARTIGO_NATURE_SP/Tree/Gene_Excisor"

SKIP_REALIGN = False
MAFFT_STRATEGY = "--auto"

RUN_IQTREE = True
IQTREE_THREADS = 30
IQTREE_ARGS = [
    "-m", "MFP",
    "-bb", "1000",
    "-alrt", "1000",
    "-nt", str(IQTREE_THREADS)
]

LOG_LEVEL = "INFO"

# ==============================================================================

import os
import sys
import subprocess
import logging
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# ==============================================================================
# LOGGING
# ==============================================================================

logging.basicConfig(
    level=getattr(logging, LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
)

log = logging.getLogger(__name__)

# ==============================================================================
# SYSTEM CHECKS
# ==============================================================================

def check_mafft():
    if subprocess.run(["which", "mafft"], capture_output=True).returncode != 0:
        log.error("MAFFT not found in PATH.")
        sys.exit(1)


def check_iqtree():
    if subprocess.run(["which", "iqtree3"], capture_output=True).returncode != 0:
        log.error("IQ-TREE 3 not found in PATH.")
        sys.exit(1)

# ==============================================================================
# UTILITIES
# ==============================================================================

def sanitize_name(name: str) -> str:
    for ch in (" ", "/", "\\", ":", "*", "?", '"', "'"):
        name = name.replace(ch, "_")
    return name.strip()

# ==============================================================================
# GENBANK PARSING
# ==============================================================================

def parse_genes(gbk_path: str):
    records = list(SeqIO.parse(gbk_path, "genbank"))
    if not records:
        log.error("No GenBank records found.")
        sys.exit(1)

    genes = []

    for feat in records[0].features:
        if feat.type not in ("CDS", "gene"):
            continue

        q = feat.qualifiers
        name = (
            q.get("locus_tag", [None])[0]
            or q.get("gene", [None])[0]
            or q.get("product", [None])[0]
            or f"feature_{feat.location.start}_{feat.location.end}"
        )

        genes.append({
            "name": sanitize_name(name),
            "start": int(feat.location.start),
            "end": int(feat.location.end),
            "strand": feat.location.strand or 1,
            "_type": feat.type,
        })

    # Deduplicate (CDS priority)
    seen = {}
    for g in sorted(genes, key=lambda x: (x["start"], x["end"])):
        key = (g["start"], g["end"])
        if key not in seen or g["_type"] == "CDS":
            seen[key] = g

    genes = sorted(seen.values(), key=lambda x: x["start"])

    log.info(f"{len(genes)} genes parsed from GenBank")
    return genes


def get_reference_from_gbk(gbk_path: str):
    record = next(SeqIO.parse(gbk_path, "genbank"))
    return SeqRecord(
        record.seq,
        id=record.id or "GBK_REF",
        description="GenBank reference"
    )

# ==============================================================================
# FASTA
# ==============================================================================

def parse_fasta(path: str):
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        log.error("Empty FASTA file.")
        sys.exit(1)
    return records

# ==============================================================================
# MAFFT REFERENCE PROJECTION
# ==============================================================================

def add_reference_to_alignment(aln_path: str, ref: SeqRecord, tmp: Path):
    ref_fa = tmp / "ref.fasta"
    out_fa = tmp / "with_ref.fasta"

    SeqIO.write([ref], ref_fa, "fasta")

    cmd = [
        "mafft",
        "--add", str(ref_fa),
        "--keeplength",
        "--quiet",
        aln_path
    ]

    log.info("Adding GenBank reference to alignment via MAFFT...")

    with open(out_fa, "w") as fh:
        result = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        log.error(result.stderr)
        sys.exit(1)

    records = list(SeqIO.parse(out_fa, "fasta"))

    # move reference to first position
    records.insert(0, records.pop(-1))

    return records

# ==============================================================================
# COORDINATE MAPPING
# ==============================================================================

def build_nt_to_col_map(seq: str):
    return [i for i, b in enumerate(seq) if b != "-"]


def extract_columns(records, nt_to_col, start, end):
    ref_len = len(nt_to_col)

    start = max(0, min(start, ref_len))
    end   = max(0, min(end, ref_len))

    if start >= end:
        return []

    col_start = nt_to_col[start]
    col_end   = nt_to_col[end - 1] + 1

    return [
        SeqRecord(
            Seq(str(r.seq)[col_start:col_end]),
            id=r.id,
            description=r.description
        )
        for r in records
    ]

# ==============================================================================
# IO
# ==============================================================================

def write_fasta(records, path):
    SeqIO.write(records, path, "fasta")


def run_mafft(in_fa, out_fa):
    cmd = ["mafft", MAFFT_STRATEGY, "--quiet", str(in_fa)]

    with open(out_fa, "w") as fh:
        result = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        log.warning(result.stderr)

# ==============================================================================
# IQ-TREE
# ==============================================================================

def run_iqtree(alignment: Path, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)

    prefix = outdir / "iqtree"

    cmd = [
        "iqtree3",
        "-s", str(alignment),
        "-pre", str(prefix),
    ] + IQTREE_ARGS

    log.info(f"Running IQ-TREE on {alignment.name}")

    result = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        log.warning(result.stderr)

# ==============================================================================
# MAIN
# ==============================================================================

def main():

    if not os.path.isfile(GBK_FILE) or not os.path.isfile(FASTA_FILE):
        log.error("Missing input files.")
        sys.exit(1)

    check_mafft()
    if RUN_IQTREE:
        check_iqtree()

    log.info("Parsing GenBank and FASTA...")

    genes = parse_genes(GBK_FILE)
    ref = get_reference_from_gbk(GBK_FILE)

    tmp = Path("_tmp_ref")
    tmp.mkdir(exist_ok=True)

    records = add_reference_to_alignment(FASTA_FILE, ref, tmp)

    nt_to_col = build_nt_to_col_map(str(records[0].seq))

    root = Path(OUTPUT_DIR)
    root.mkdir(exist_ok=True)

    skipped = 0

    log.info("Processing genes...")

    for g in genes:
        name, start, end = g["name"], g["start"], g["end"]

        log.info(f"{name} [{start}:{end}]")

        sliced = extract_columns(records, nt_to_col, start, end)

        if not sliced:
            skipped += 1
            log.warning("Skipped (out of range)")
            continue

        sliced = sliced[1:]  # remove reference

        gene_dir = root / name
        gene_dir.mkdir(exist_ok=True)

        raw = gene_dir / f"{name}_raw.fasta"
        aln = gene_dir / f"{name}_aligned.fasta"

        write_fasta(sliced, raw)

        if SKIP_REALIGN:
            aln = raw
        else:
            run_mafft(raw, aln)

        if RUN_IQTREE:
            run_iqtree(aln, gene_dir / "iqtree")

    log.info(f"Done: {len(genes)-skipped}/{len(genes)} genes processed")


if __name__ == "__main__":
    main()