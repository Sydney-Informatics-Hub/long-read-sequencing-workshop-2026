#!/usr/bin/env python3
"""
contig_diagram.py

Draw a diagram of one or more contigs (each linear or circular) with
regions of interest overlaid as strand-oriented arrows. Coverage controls
arrow colour (a viridis gradient); percentage identity controls arrow
opacity (faint = low identity, solid = high/100% identity).

Why two libraries?
-------------------
pyGenomeViz only draws LINEAR tracks. There is no native "circular contig"
mode in pyGenomeViz itself. For circular replicons (plasmids, circularised
chromosomes, etc.) this script uses pyCirclize (the circular-plot sibling
package from the same author) so that circularity is actually represented
in the figure rather than just labelled. Linear contigs are drawn with
pyGenomeViz. Both use the same coverage colour scale so panels are
visually comparable, and are composited into one figure with a shared
colourbar.

Install:
    pip install pygenomeviz pycirclize matplotlib

Usage:
    Edit CONTIGS and REGIONS below (or load them from a CSV/TSV -- see
    load_contigs_csv / load_regions_csv), then run:
        python contig_diagram.py
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Literal
from sys import argv

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.image import imread
import numpy as np

from pygenomeviz import GenomeViz
from pycirclize import Circos


# --------------------------------------------------------------------------
# Data model
# --------------------------------------------------------------------------

@dataclass
class Contig:
    name: str
    length: int
    topology: Literal["circular", "linear"]


@dataclass
class Region:
    contig: str
    name: str
    start: int
    end: int
    strand: Literal[1, -1]
    coverage: float
    pident: float = 100.0  # percentage identity, 0-100

# --------------------------------------------------------------------------
# Optional: load from CSV/TSV instead of hardcoding
# --------------------------------------------------------------------------

def load_contigs_csv(path: str) -> list[Contig]:
    """Expected columns: name, length, topology (circular/linear)."""
    contigs = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            contigs.append(
                Contig(
                    name=row["name"],
                    length=int(row["length"]),
                    topology=row["topology"].strip().lower(),
                )
            )
    return contigs


def load_regions_csv(path: str) -> list[Region]:
    """Expected columns: contig, name, start, end, strand, coverage, pident.

    `pident` (percentage identity, 0-100) is optional in the file; if the
    column is missing or a row's value is blank, it defaults to 100 (fully
    opaque)."""
    regions = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            pident_raw = row.get("pident", "")
            pident = float(pident_raw) if pident_raw not in (None, "") else 100.0
            regions.append(
                Region(
                    contig=row["contig"],
                    name=row["name"],
                    start=int(row["start"]),
                    end=int(row["end"]),
                    strand=int(row["strand"]),
                    coverage=float(row["coverage"]),
                    pident=pident,
                )
            )
    return regions


# --------------------------------------------------------------------------
# Colour mapping
# --------------------------------------------------------------------------

def make_coverage_cmap(regions: list[Region], cmap_name: str = "viridis"):
    """Build a shared colormap + normaliser across ALL regions so every
    panel (linear or circular) uses the same coverage colour scale."""
    coverages = [r.coverage for r in regions]
    cmap = plt.get_cmap(cmap_name)
    norm = mcolors.Normalize(vmin=min(coverages), vmax=max(coverages))
    return cmap, norm


def coverage_to_hex(cov: float, cmap, norm) -> str:
    return mcolors.to_hex(cmap(norm(cov)))


# Percentage-identity -> opacity mapping. Identity is treated as an
# absolute 0-100 scale (not rescaled to the data's min/max) so that, e.g.,
# 100% identity always means fully opaque and figures stay comparable
# across different datasets. MIN_ALPHA keeps even low-identity regions
# faintly visible rather than disappearing entirely.
MIN_ALPHA = 0.15
MAX_ALPHA = 1.0
PIDENT_VMIN = 0.0
PIDENT_VMAX = 100.0


def pident_to_alpha(pident: float) -> float:
    pident = max(PIDENT_VMIN, min(PIDENT_VMAX, pident))
    frac = (pident - PIDENT_VMIN) / (PIDENT_VMAX - PIDENT_VMIN)
    return MIN_ALPHA + frac * (MAX_ALPHA - MIN_ALPHA)


def region_facecolor(region: Region, cmap, norm) -> tuple[float, float, float, float]:
    """RGBA colour for a region: hue/lightness from coverage, opacity from
    percentage identity. Returned as an explicit RGBA tuple (rather than
    using the patch-level `alpha` kwarg) so that opacity affects only the
    fill and not the edge, keeping black outlines crisp at any identity."""
    r, g, b, _ = cmap(norm(region.coverage))
    return (r, g, b, pident_to_alpha(region.pident))


# --------------------------------------------------------------------------
# Linear contig rendering (pyGenomeViz)
# --------------------------------------------------------------------------

def render_linear_contig(contig: Contig, regions: list[Region], cmap, norm,
                          outpath: Path) -> Path:
    gv = GenomeViz(fig_track_height=1.2, feature_track_ratio=0.6)
    gv.set_scale_bar()
    track = gv.add_feature_track(contig.name, contig.length)
    track.add_sublabel(f"{contig.length:,} bp (linear)")

    for r in regions:
        color = region_facecolor(r, cmap, norm)
        track.add_feature(
            r.start, r.end, r.strand,
            plotstyle="bigarrow",
            fc=color, ec="black", lw=0.4,
            label=r.name,
            text_kws=dict(size=8, rotation=30, vpos="top", hpos="left"),
        )

    fig = gv.plotfig()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return outpath


# --------------------------------------------------------------------------
# Circular contig rendering (pyCirclize)
# --------------------------------------------------------------------------

def render_circular_contig(contig: Contig, regions: list[Region], cmap, norm,
                            outpath: Path) -> Path:
    circos = Circos(sectors={contig.name: contig.length}, space=0)
    sector = circos.sectors[0]

    circos.text(f"{contig.name}\n{contig.length:,} bp (circular)",
                r=20, size=11)

    # Outer axis ring with tick marks
    axis_track = sector.add_track((98, 100))
    axis_track.axis(fc="lightgrey", ec="none")
    major = max(1000, int(10 ** math.floor(math.log10(contig.length / 4))))
    axis_track.xticks_by_interval(
        major, label_formatter=lambda v: f"{v / 1000:.1f} kb"
    )

    # Separate tracks for forward (outer) and reverse (inner) strand arrows
    fwd_track = sector.add_track((90, 97))
    rev_track = sector.add_track((82, 89))

    for r in regions:
        color = region_facecolor(r, cmap, norm)
        target = fwd_track if r.strand == 1 else rev_track
        # Arrowhead points from `start` towards `end`; for reverse-strand
        # features we swap the coordinates so the arrow points "backwards".
        arrow_start, arrow_end = (r.start, r.end) if r.strand == 1 else (r.end, r.start)
        target.arrow(
            arrow_start, arrow_end,
            fc=color, ec="black", lw=0.4,
        )
        # Label near the midpoint of the feature
        mid = (r.start + r.end) / 2
        target.text(r.name, x=mid, size=9, orientation="vertical", adjust_rotation=True)

    fig = circos.plotfig()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return outpath


# --------------------------------------------------------------------------
# Compose everything into one figure with a shared colourbar
# --------------------------------------------------------------------------

def build_diagram(contigs: list[Contig], regions: list[Region],
                   outfile: str = "contig_diagram.png",
                   cmap_name: str = "viridis") -> Path:
    cmap, norm = make_coverage_cmap(regions, cmap_name)
    regions_by_contig = {c.name: [] for c in contigs}
    for r in regions:
        if r.contig not in regions_by_contig:
            raise ValueError(f"Region '{r.name}' references unknown contig '{r.contig}'")
        regions_by_contig[r.contig].append(r)

    tmp_dir = Path("panel_tmp")
    tmp_dir.mkdir(exist_ok=True)
    panel_paths = []
    for contig in contigs:
        creg = regions_by_contig.get(contig.name, [])
        panel_path = tmp_dir / f"{contig.name}.png"
        if contig.topology == "linear":
            render_linear_contig(contig, creg, cmap, norm, panel_path)
        elif contig.topology == "circular":
            render_circular_contig(contig, creg, cmap, norm, panel_path)
        else:
            raise ValueError(f"Unknown topology '{contig.topology}' for contig '{contig.name}'")
        panel_paths.append(panel_path)

    # Composite panels vertically, plus shared legends at the bottom:
    # a colourbar for coverage and a greyscale opacity strip for identity.
    n = len(panel_paths)
    fig, axes = plt.subplots(
        n + 2, 1, figsize=(9, 3.2 * n + 1.1),
        gridspec_kw={"height_ratios": [1] * n + [0.06, 0.06]},
    )

    for ax, panel_path in zip(axes[:n], panel_paths):
        ax.imshow(imread(panel_path))
        ax.axis("off")

    cov_cbar_ax = axes[n]
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, cax=cov_cbar_ax, orientation="horizontal", label="Coverage")

    # Opacity legend: black fading from MIN_ALPHA to MAX_ALPHA across the
    # 0-100% identity range, matching how region_facecolor() sets alpha.
    identity_ax = axes[n + 1]
    n_steps = 256
    gradient = np.zeros((1, n_steps, 4))
    gradient[0, :, 3] = np.linspace(MIN_ALPHA, MAX_ALPHA, n_steps)  # alpha channel
    identity_ax.imshow(
        gradient, aspect="auto",
        extent=(PIDENT_VMIN, PIDENT_VMAX, 0, 1),
    )
    identity_ax.set_yticks([])
    identity_ax.set_xlim(PIDENT_VMIN, PIDENT_VMAX)
    identity_ax.set_xlabel("Percentage identity (opacity)")

    fig.tight_layout()
    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved combined diagram to {outfile}")
    return Path(outfile)


if __name__ == "__main__":
    assert len(argv) == 4, 'Usage: python3 contig_diagram.py <CONTIGS>.csv <REGIONS>.csv <OUTPUT>.png.'
    contigs = load_contigs_csv(argv[1])
    regions = load_regions_csv(argv[2])
    output = argv[3]
    build_diagram(contigs, regions, outfile=output)
