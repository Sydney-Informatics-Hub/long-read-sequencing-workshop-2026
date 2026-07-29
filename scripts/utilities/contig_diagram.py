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
import sys

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.image import imread
import numpy as np

from pygenomeviz import GenomeViz
from pycirclize import Circos


# --------------------------------------------------------------------------
# Global sizing -- tune these to make tracks/text bigger or smaller
# everywhere at once.
# --------------------------------------------------------------------------

# pyGenomeViz (linear contigs)
LINEAR_FIG_TRACK_HEIGHT = 2.2     # overall height of each linear track
LINEAR_FEATURE_TRACK_RATIO = 0.7  # fraction of track height used by the arrow itself
LINEAR_TRACK_LABELSIZE = 24       # contig name, e.g. "chromosome" (pyGenomeViz default: 20)
LINEAR_SUBLABEL_SIZE = 24         # "42,000 bp (linear)" line
LINEAR_FEATURE_LABEL_SIZE = 24    # gene/region name labels
LINEAR_SCALEBAR_LABELSIZE = 20    # scale bar text
LINEAR_BACKBONE_LW = 5.0          # thickness of the horizontal contig line

# pyCirclize (circular contigs)
CIRCULAR_FIGSIZE = (10, 10)       # pyCirclize default: (8, 8)
CIRCULAR_TITLE_SIZE = 24          # contig name + length, e.g. "plasmid_1\n9,200 bp"
CIRCULAR_TICK_LABEL_SIZE = 20     # kb tick labels around the ring
CIRCULAR_FEATURE_LABEL_SIZE = 24  # gene/region name labels
CIRCULAR_FWD_TRACK_R = (83, 95)   # forward-strand track's full radius range, wider = thicker arrows
CIRCULAR_REV_TRACK_R = (68, 80)   # reverse-strand track's full radius range
CIRCULAR_AXIS_R = (96, 100)       # grey backbone ring radius range, wider = thicker
# Each strand track above is split into an outer zone for the arrow itself
# and an inner zone for its label, so labels sit inward of the bars rather
# than on top of them. Fractions are of that track's own radius range.
CIRCULAR_LABEL_ZONE_FRAC = 0.4    # inner 40% is reserved for the label
CIRCULAR_MAX_TICKS = 12           # cap on axis ticks so labels don't crowd large contigs

# Composite figure (final combined PNG)
PANEL_HEIGHT_PER_CONTIG = 4.5     # inches of vertical space per contig panel (was 3.2)
COMPOSITE_DPI = 200               # was 200
LEGEND_LABEL_FONTSIZE = 13        # "Coverage" / "Percentage identity (opacity)" labels
LEGEND_TICK_FONTSIZE = 11         # numbers along the legend bars


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


def nice_tick_interval(length: int, max_ticks: int) -> int:
    """Pick a 'round' tick interval (1/2/5 x a power of ten) so that the
    number of ticks around a contig of the given length stays at or below
    max_ticks. This keeps kb labels from overlapping on large contigs."""
    raw_interval = length / max_ticks
    magnitude = 10 ** math.floor(math.log10(raw_interval))
    for step in (1, 2, 5, 10):
        candidate = step * magnitude
        if raw_interval <= candidate:
            return max(1, int(candidate))
    return max(1, int(10 * magnitude))


# --------------------------------------------------------------------------
# Linear contig rendering (pyGenomeViz)
# --------------------------------------------------------------------------

def render_linear_contig(contig: Contig, regions: list[Region], cmap, norm,
                          outpath: Path) -> Path:
    gv = GenomeViz(
        fig_track_height=LINEAR_FIG_TRACK_HEIGHT,
        feature_track_ratio=LINEAR_FEATURE_TRACK_RATIO,
    )
    gv.set_scale_bar(labelsize=LINEAR_SCALEBAR_LABELSIZE)
    track = gv.add_feature_track(
        contig.name, contig.length, labelsize=LINEAR_TRACK_LABELSIZE,
        line_kws=dict(lw=LINEAR_BACKBONE_LW, color="grey"),
    )
    track.add_sublabel(
        f"{contig.length:,} bp (linear)", size=LINEAR_SUBLABEL_SIZE
    )

    for r in regions:
        color = region_facecolor(r, cmap, norm)
        track.add_feature(
            r.start, r.end, r.strand,
            plotstyle="bigarrow",
            fc=color, ec="black", lw=0.4,
            label=r.name,
            text_kws=dict(size=LINEAR_FEATURE_LABEL_SIZE, rotation=30,
                           vpos="top", hpos="left"),
        )

    fig = gv.plotfig()
    fig.savefig(outpath, dpi=COMPOSITE_DPI, bbox_inches="tight")
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
                r=20, size=CIRCULAR_TITLE_SIZE)

    # Outer axis ring with tick marks. Interval is chosen so the number of
    # ticks stays at or below CIRCULAR_MAX_TICKS, regardless of contig size.
    axis_track = sector.add_track(CIRCULAR_AXIS_R)
    axis_track.axis(fc="lightgrey", ec="none")
    interval = nice_tick_interval(contig.length, CIRCULAR_MAX_TICKS)
    axis_track.xticks_by_interval(
        interval, label_formatter=lambda v: f"{v / 1000:.1f} kb",
        label_size=CIRCULAR_TICK_LABEL_SIZE,
    )

    # Separate tracks for forward (outer) and reverse (inner) strand
    # features. Each track is split radially into an outer "arrow zone"
    # (where the strand arrow is drawn) and an inner "label zone" (where
    # the feature name sits), so labels don't overlap the bars.
    fwd_track = sector.add_track(CIRCULAR_FWD_TRACK_R)
    rev_track = sector.add_track(CIRCULAR_REV_TRACK_R)

    def zone_radii(track) -> tuple[tuple[float, float], float]:
        r0, r1 = track.r_lim
        span = r1 - r0
        arrow_r_lim = (r0 + span * CIRCULAR_LABEL_ZONE_FRAC, r1)
        label_r = r0 + span * CIRCULAR_LABEL_ZONE_FRAC * 0.1
        return arrow_r_lim, label_r

    fwd_arrow_r, fwd_label_r = zone_radii(fwd_track)
    rev_arrow_r, rev_label_r = zone_radii(rev_track)

    for r in regions:
        color = region_facecolor(r, cmap, norm)
        if r.strand == 1:
            target, arrow_r_lim, label_r = fwd_track, fwd_arrow_r, fwd_label_r
        else:
            target, arrow_r_lim, label_r = rev_track, rev_arrow_r, rev_label_r

        # Arrowhead points from `start` towards `end`; for reverse-strand
        # features we swap the coordinates so the arrow points "backwards".
        arrow_start, arrow_end = (r.start, r.end) if r.strand == 1 else (r.end, r.start)
        target.arrow(
            arrow_start, arrow_end,
            r_lim=arrow_r_lim,
            fc=color, ec="black", lw=0.4,
        )
        # Label offset inward from the arrow, oriented tangentially so it
        # runs along the length of the feature bar rather than radially.
        mid = (r.start + r.end) / 2
        target.text(r.name, x=mid, r=label_r,
                     size=CIRCULAR_FEATURE_LABEL_SIZE,
                     orientation="horizontal", adjust_rotation=True)

    fig = circos.plotfig(figsize=CIRCULAR_FIGSIZE, dpi=COMPOSITE_DPI)
    fig.savefig(outpath, dpi=COMPOSITE_DPI, bbox_inches="tight")
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
        n + 2, 1, figsize=(9, PANEL_HEIGHT_PER_CONTIG * n + 1.4),
        gridspec_kw={"height_ratios": [1] * n + [0.07, 0.07]},
    )

    for ax, panel_path in zip(axes[:n], panel_paths):
        ax.imshow(imread(panel_path))
        ax.axis("off")

    cov_cbar_ax = axes[n]
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cov_cbar = fig.colorbar(sm, cax=cov_cbar_ax, orientation="horizontal")
    cov_cbar.set_label("Coverage", fontsize=LEGEND_LABEL_FONTSIZE)
    cov_cbar.ax.tick_params(labelsize=LEGEND_TICK_FONTSIZE)

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
    identity_ax.set_xlabel("Percentage identity (opacity)", fontsize=LEGEND_LABEL_FONTSIZE)
    identity_ax.tick_params(labelsize=LEGEND_TICK_FONTSIZE)

    fig.tight_layout()
    fig.savefig(outfile, dpi=COMPOSITE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved combined diagram to {outfile}")
    return Path(outfile)


if __name__ == "__main__":
    assert len(sys.argv) == 4, 'Usage: python3 contig_diagram.py <CONTIGS>.csv <REGIONS>.csv <OUTPUT>.png.'
    contigs = load_contigs_csv(sys.argv[1])
    regions = load_regions_csv(sys.argv[2])
    output = sys.argv[3]
    build_diagram(contigs, regions, outfile=output)
