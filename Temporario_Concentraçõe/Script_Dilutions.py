#!/usr/bin/env python3
"""
Dilution calculator: water needed to bring samples to a target concentration.

Uses C1*V1 = C2*V2.
Edit the CONFIG block below and run (VS Code: F5 / "Run Python File").
"""

# =========================== CONFIG ===========================

TARGET_CONC = 20.0        # ng/µl, desired final concentration

# MODE decides which volumes you fix:
#   "sample_volume"    -> fixed sample volume (SAMPLE_VOLUME_UL);
#                         returns water and resulting final volume.
#   "final_volume"     -> fixed final volume (FINAL_VOLUME_UL);
#                         returns sample volume and water.
#   "max_sample_final" -> desired final volume (FINAL_VOLUME_UL) with a cap on
#                         sample usage (MAX_SAMPLE_VOLUME_UL). If the required
#                         sample exceeds the cap, the cap is used and the final
#                         volume is reduced accordingly (flagged in the table).
MODE = "max_sample_final"

SAMPLE_VOLUME_UL = 5.0        # µl, MODE "sample_volume"
FINAL_VOLUME_UL = 10.0        # µl, MODE "final_volume" / "max_sample_final"
MAX_SAMPLE_VOLUME_UL = 10.0   # µl, MODE "max_sample_final"

# Initial concentrations in ng/µl, one entry per sample
SAMPLES = {
    "S01": 37.7,
    "S02": 36.5,
    "S03": 24.2,
    "S04": 20.9,    # below target: cannot be diluted
    "S05": 27.7,
    "S06": 28.3,
    "S07": 41.0,
    "S08": 15.6,
    "S09": 1.64,    # below target: cannot be diluted
    "S10": 38.5,
    "S11": 28.6,
    "S12": 12.8,
    "S13": 29.8,
    "S14": 14.8,    # below target: cannot be diluted
    "S15": 0.532,
    "S16": 32.0,
    "S17": 36.3,
    "S18": 34.5,
    "S19": 24.1,    # below target: cannot be diluted
    "S20": 38.8,
    "S21": 31.1,
    "S22": 42.2,
    "S23": 20.9,
    "S24": 33.1,    # below target: cannot be diluted     
}

# Optional: save results as CSV (set to None to disable)
OUTPUT_CSV = "dilutions_16S.csv"         # e.g. "dilutions.csv"

# ===============================================================


def calc_dilution(conc, target):
    """Return (sample_ul, water_ul, final_ul, note) for one sample."""
    if conc <= 0:
        return None, None, None, "invalid concentration"
    if conc < target:
        return None, None, None, "below target, no dilution (concentrate or use as is)"

    note = "already at target" if conc == target else ""

    if MODE == "sample_volume":
        sample_ul = SAMPLE_VOLUME_UL
        final_ul = conc * sample_ul / target

    elif MODE == "final_volume":
        final_ul = FINAL_VOLUME_UL
        sample_ul = target * final_ul / conc

    elif MODE == "max_sample_final":
        final_ul = FINAL_VOLUME_UL
        sample_ul = target * final_ul / conc
        if sample_ul > MAX_SAMPLE_VOLUME_UL:
            sample_ul = MAX_SAMPLE_VOLUME_UL
            final_ul = conc * sample_ul / target
            note = f"sample capped; final volume reduced from {FINAL_VOLUME_UL:g} µl"
    else:
        raise ValueError(f"Unknown MODE: {MODE!r}")

    return sample_ul, final_ul - sample_ul, final_ul, note


def main():
    rows = []
    for name, conc in SAMPLES.items():
        rows.append((name, conc, *calc_dilution(conc, TARGET_CONC)))

    fmt = lambda x: "-" if x is None else f"{x:.2f}"

    print(f"Target: {TARGET_CONC:g} ng/µl | Mode: {MODE}\n")
    header = f"{'Sample':<10}{'Conc (ng/µl)':>14}{'Sample (µl)':>13}{'Water (µl)':>12}{'Final (µl)':>12}  Note"
    print(header)
    print("-" * len(header))
    for name, conc, s, w, f, note in rows:
        print(f"{name:<10}{conc:>14.2f}{fmt(s):>13}{fmt(w):>12}{fmt(f):>12}  {note}")

    total_water = sum(r[3] for r in rows if r[3] is not None)
    print(f"\nTotal water: {total_water:.2f} µl")

    if OUTPUT_CSV:
        import csv
        with open(OUTPUT_CSV, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["sample", "conc_ng_ul", "sample_ul", "water_ul", "final_ul", "note"])
            for name, conc, s, w, f, note in rows:
                writer.writerow([name, conc, fmt(s), fmt(w), fmt(f), note])
        print(f"Saved: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()