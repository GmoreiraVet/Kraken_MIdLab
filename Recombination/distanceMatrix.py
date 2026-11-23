#!/usr/bin/env python3
"""Clustered distance heatmaps from IQ-TREE .mldist (PHYLIP square) files.

One PNG per input file. Each heatmap has a square plot area, and the figure
size is computed from the number of taxa so the tiles stay square and the
whitespace is minimal. Taxon names are reduced to the accession number.
"""

import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from plotnine import *
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

# ───────────────────────────── CONFIG ─────────────────────────────
# Point directly to the .mldist file(s). Each entry is either a path
# or a (panel tag, path) tuple. One entry gives a single heatmap.
PANELS = [
    "/home/viroicbas2023/Documents/Gmoreira/Carraças/Carraças_Metagenomica_Gui/VOLZSHKOESEPARADO/Filogenia/L_Segment/L_Segment_Database_ALIGNED_WithPTSequence.fasta.mldist",
    "/home/viroicbas2023/Documents/Gmoreira/Carraças/Carraças_Metagenomica_Gui/VOLZSHKOESEPARADO/Filogenia/M_Segment/M_Segment_Database_FInal_ConsensusPT.fasta.mldist",
    # ("B", "/path/to/second_file.mldist"),
]
OUT_DIR = Path(".")                          # where the PNGs are written
OUT_PATTERN = "{stem}_heatmap.png"           # {stem} = input file name, {tag} = panel tag
DPI = 400

CLEAN_NAMES = True      # reduce taxon names to the accession (set False to keep full names)
TILE_INCHES = 0.8       # side of one tile; the plot area is n_taxa * TILE_INCHES square

LINKAGE_METHOD = "single"  # "single", "complete", "average", "weighted", "centroid", "median", "ward"
FILL_LIMITS = (0, 1)  # distances outside this range are drawn as NA (grey)
FILL_LOW, FILL_HIGH = "#3AA288", "#13434E"
TILE_COLOR = TEXT_COLOR = "#F0E5D2"
LABEL_DECIMALS = 2
LABEL_SIZE = 20
LEGEND_TITLE_SIZE = 24   # "Distance"
LEGEND_TEXT_SIZE = 20    # tick numbers on the bar
LEGEND_BAR_WIDTH = 40    # bar thickness
LEGEND_BAR_HEIGHT = 500  # bar length
# ──────────────────────────────────────────────────────────────────

# GenBank-style (e.g. PP473488.1, OQ849227.1) or RefSeq-style (e.g. NC_012345.1) accession.version
ACCESSION_RE = re.compile(r"(?<![A-Za-z0-9])(?:[A-Z]{1,2}\d{5,8}|[A-Z]{2}_\d{6,9})\.\d+")


def clean_name(name: str) -> str:
    """Keep only the accession.version of a taxon name.

    Names without a recognisable accession (e.g. consensus or in-house
    sequences) fall back to the text before the first '/' or '|'.
    """
    m = ACCESSION_RE.search(name)
    if m:
        return m.group(0)
    return re.split(r"[/|]", name, maxsplit=1)[0].strip("_ ") or name


def clean_names(names: list[str]) -> list[str]:
    """Clean all names and make them unique (suffix _2, _3, ... on collisions)."""
    cleaned = [clean_name(n) for n in names]
    totals = Counter(cleaned)
    seen: Counter = Counter()
    unique = []
    for c in cleaned:
        seen[c] += 1
        unique.append(c if totals[c] == 1 and seen[c] == 1 else f"{c}_{seen[c]}")
    dups = [c for c, k in totals.items() if k > 1]
    if dups:
        print(f"  warning: {len(dups)} accession(s) occur more than once, suffixed: "
              f"{', '.join(dups[:5])}{'...' if len(dups) > 5 else ''}")
    return unique


def read_mldist(path: Path) -> pd.DataFrame:
    """Parse a PHYLIP square distance matrix (first line = number of taxa)."""
    lines = [l for l in path.read_text().splitlines() if l.strip()]
    n = int(lines[0])
    if len(lines) - 1 != n:
        raise ValueError(f"{path.name}: header says {n} taxa, found {len(lines) - 1} rows")

    names, rows = [], []
    for line in lines[1:]:
        name, *values = line.split()
        names.append(name)
        rows.append(np.asarray(values, dtype=float))

    if len(set(names)) != n:
        raise ValueError(f"{path.name}: duplicate taxon names")
    if any(len(r) != n for r in rows):
        raise ValueError(f"{path.name}: matrix is not square")

    if CLEAN_NAMES:
        names = clean_names(names)

    return pd.DataFrame(np.vstack(rows), index=names, columns=names)


def to_plot_frame(m: pd.DataFrame) -> pd.DataFrame:
    """Cluster taxa, then return a long frame with ordered categorical axes."""
    sym = (m + m.T) / 2                      # absorb tiny numerical asymmetry
    Z = linkage(squareform(sym.values, checks=False), method=LINKAGE_METHOD)
    order = sym.index[leaves_list(Z)].tolist()

    long = (
        sym.loc[order, order]
        .rename_axis(index="taxon_y", columns="taxon_x")
        .stack()
        .rename("distance")
        .reset_index()
    )
    long["taxon_x"] = pd.Categorical(long["taxon_x"], categories=order)
    long["taxon_y"] = pd.Categorical(long["taxon_y"], categories=order[::-1])
    long["label"] = long["distance"].map(f"{{:.{LABEL_DECIMALS}f}}".format)
    return long


def figure_size(names) -> tuple[float, float]:
    """Figure size (in) = square plot area + room for axis labels and legend."""
    n = len(names)
    panel = n * TILE_INCHES
    label_in = max(len(s) for s in names) * LABEL_SIZE * 0.6 / 72   # longest taxon name
    legend_in = (LEGEND_BAR_WIDTH + LEGEND_TEXT_SIZE * 4) / 72 + 1.0
    width = panel + label_in + legend_in + 0.5   # + y labels, legend, margin
    height = panel + label_in + 0.5              # + rotated x labels, margin
    return width, height


def make_heatmap(long: pd.DataFrame, tag: str, size: tuple[float, float]):
    return (
        ggplot(long, aes("taxon_x", "taxon_y", fill="distance"))
        + geom_tile(color=TILE_COLOR, size=0.3)
        + geom_text(aes(label="label"), color=TEXT_COLOR, size=LABEL_SIZE)
        + scale_fill_gradient(low=FILL_LOW, high=FILL_HIGH, limits=FILL_LIMITS, name="Distance")
        + labs(tag=tag)
        + coord_fixed(ratio=1, expand=False)        # square tiles -> square plot area
        + theme_minimal()
        + theme(
            axis_text_x=element_text(rotation=90, ha="center", size=LABEL_SIZE),
            axis_text_y=element_text(size=LABEL_SIZE),
            axis_title=element_blank(),
            panel_grid=element_blank(),
            legend_title=element_text(size=LEGEND_TITLE_SIZE),
            legend_text=element_text(size=LEGEND_TEXT_SIZE),
            legend_key_width=LEGEND_BAR_WIDTH,
            legend_key_height=LEGEND_BAR_HEIGHT,
            figure_size=size,
        )
    )


def normalise_panels(panels):
    """Accept paths or (tag, path) tuples; return a list of (tag, Path)."""
    out = []
    for entry in panels:
        tag, path = ("", entry) if isinstance(entry, (str, Path)) else entry
        out.append((tag, Path(path).expanduser()))
    return out


def main() -> None:
    panels = normalise_panels(PANELS)
    if not panels:
        raise SystemExit("PANELS is empty: add at least one .mldist file.")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for tag, path in panels:
        if not path.is_file():
            raise FileNotFoundError(path)
        m = read_mldist(path)
        print(f"{tag or path.stem}: {path.name} | taxa: {len(m)} | "
              f"range: {m.values.min():.4f}–{m.values.max():.4f}")

        size = figure_size(list(m.index))
        p = make_heatmap(to_plot_frame(m), tag, size)

        out = OUT_DIR / OUT_PATTERN.format(stem=path.name.removesuffix(".mldist"), tag=tag)
        p.save(out, dpi=DPI, limitsize=False)
        print(f"  saved {out} ({size[0]:.1f} x {size[1]:.1f} in)")


if __name__ == "__main__":
    main()