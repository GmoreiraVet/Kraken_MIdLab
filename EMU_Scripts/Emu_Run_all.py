#!/usr/bin/env python3
"""Run Emu on every FASTQ/FASTA file in a folder, then merge results.

Edit the CONFIG block below and press Run in VS Code (no arguments needed).
"""

import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# =============================== CONFIG ===============================
INPUT_DIR = Path("/home/viroicbas2023/Documents/Gmoreira/AbelhasSoraia/Plamid16S_Outubro26/Reads/QV3LQH_fastq")
OUTPUT_DIR = Path("/home/viroicbas2023/Documents/Gmoreira/AbelhasSoraia/emu_results")
DB_DIR = Path("/home/viroicbas2023/Documents/DATABASES/EMU/emu")        # must contain species_taxid.fasta + taxonomy.tsv

# Path to the emu executable. "emu" works if the conda env is active in the
# VS Code terminal/interpreter; otherwise give the full path, e.g.
# "/home/you/miniconda3/envs/emu/bin/emu"
EMU_CMD = "emu"

EMU_TYPE = "lr:hq"    # map-ont | lr:hq | map-pb | sr
THREADS = 30             # threads per sample
JOBS = 1                # samples run in parallel (total CPU = JOBS x THREADS)
KEEP_COUNTS = True      # add estimated read counts to output
EXTRA_ARGS = []         # e.g. ["--min-abundance", "0.001", "--keep-read-assignments"]
SKIP_EXISTING = True    # skip samples that already have an output file
COMBINE_RANKS = ["species", "genus"]   # merged tables at these ranks; [] to disable
COMBINE_COUNTS = True   # also write a species-level counts table
# ======================================================================

READ_RE = re.compile(r"\.(fastq|fq|fasta|fa)(\.gz)?$", re.IGNORECASE)


def sample_name(path: Path) -> str:
    return READ_RE.sub("", path.name)


def run_sample(read_file: Path):
    sample = sample_name(read_file)
    result = OUTPUT_DIR / f"{sample}_rel-abundance.tsv"
    log = OUTPUT_DIR / "logs" / f"{sample}.log"

    if SKIP_EXISTING and result.exists():
        return sample, "skipped", None

    cmd = [
        EMU_CMD, "abundance", str(read_file),
        "--db", str(DB_DIR),
        "--type", EMU_TYPE,
        "--threads", str(THREADS),
        "--output-dir", str(OUTPUT_DIR),
        "--output-basename", sample,
    ]
    if KEEP_COUNTS:
        cmd.append("--keep-counts")
    cmd += EXTRA_ARGS

    with open(log, "w") as fh:
        fh.write("CMD: " + " ".join(cmd) + "\n\n")
        fh.flush()
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    return sample, ("done" if proc.returncode == 0 else "failed"), log


def combine(rank: str, counts: bool = False):
    cmd = [EMU_CMD, "combine-outputs", str(OUTPUT_DIR), rank]
    if counts:
        cmd.append("--counts")
    log = OUTPUT_DIR / "logs" / f"combine_{rank}{'_counts' if counts else ''}.log"
    with open(log, "w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        print(f"WARN: combine-outputs ({rank}) failed, see {log}", file=sys.stderr)


def main():
    if shutil.which(EMU_CMD) is None and not Path(EMU_CMD).is_file():
        sys.exit(f"ERROR: '{EMU_CMD}' not found. Activate the emu conda env or set EMU_CMD.")
    if not INPUT_DIR.is_dir():
        sys.exit(f"ERROR: INPUT_DIR not found: {INPUT_DIR}")
    for required in ("species_taxid.fasta", "taxonomy.tsv"):
        if not (DB_DIR / required).is_file():
            sys.exit(f"ERROR: {DB_DIR} is missing {required}")

    files = sorted(f for f in INPUT_DIR.iterdir() if f.is_file() and READ_RE.search(f.name))
    if not files:
        sys.exit(f"ERROR: no read files found in {INPUT_DIR}")

    (OUTPUT_DIR / "logs").mkdir(parents=True, exist_ok=True)
    print(f"Found {len(files)} files | db={DB_DIR} | type={EMU_TYPE} | "
          f"threads={THREADS} | jobs={JOBS}")

    failed = []
    with ThreadPoolExecutor(max_workers=JOBS) as pool:
        futures = [pool.submit(run_sample, f) for f in files]
        for fut in as_completed(futures):
            sample, status, log = fut.result()
            if status == "failed":
                failed.append(sample)
                print(f"[FAIL] {sample} (see {log})", file=sys.stderr)
            else:
                print(f"[{status}] {sample}")

    if any(OUTPUT_DIR.glob("*_rel-abundance.tsv")):
        for rank in COMBINE_RANKS:
            combine(rank)
        if COMBINE_COUNTS and KEEP_COUNTS:
            combine("species", counts=True)

    if failed:
        sys.exit(f"Finished with failures: {', '.join(failed)}")
    print(f"All samples finished. Results in {OUTPUT_DIR}")


if __name__ == "__main__":
    main()