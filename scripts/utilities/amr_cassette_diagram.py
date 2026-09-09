#!/usr/bin/env python3
"""
amr_cassette_diagram.py

Draw a linear gene-cassette map straight from an NCBI AMRFinderPlus TSV
report: one contig's resistance/stress hits as strand-oriented arrows,
positioned at their real coordinates, coloured by antibiotic class
(categorical), with arrow opacity showing %coverage of the reference gene
-- so a PARTIALX hit visibly fades relative to a full ALLELEX/EXACTX one.
Optionally composites a pre-rendered Bandage assembly-graph image
underneath as a second panel, to show whether the cassette-carrying
contig is the chromosome or a plasmid.

Companion to contig_diagram.py in this same folder: same pyGenomeViz
linear-track machinery, but the AMRFinderPlus TSV is read directly (no
separate CSV conversion needed), and colour encodes drug class rather
than BLAST coverage, since class is the informative axis for an AMR hit.

Install (conda/mamba):
    conda create -n amrviz -c bioconda -c conda-forge pygenomeviz
    conda activate amrviz

Or as a BioContainer (pyGenomeViz is a bioconda package, so it ships as
one -- no local install needed; confirmed present at
quay.io/biocontainers/pygenomeviz, tag 0.4.4--pyhdfd78af_0 as of writing).

IMPORTANT: bioconda/biocontainers has been stuck on pyGenomeViz 0.4.4 since
2023 -- pyGenomeViz did a breaking 1.0 API rewrite upstream (now at 1.7.0 on
PyPI) that bioconda's recipe was never updated to track. This script is
written against the 0.4.4 API on purpose, to match what the container
actually provides (`pip install pygenomeviz` outside a container will give
you the newer, incompatible 1.x API instead -- don't mix the two). Check
quay.io/repository/biocontainers/pygenomeviz?tab=tags in case that's
changed:
    docker run --rm -v "$PWD":/data -w /data \\
        quay.io/biocontainers/pygenomeviz:0.4.4--pyhdfd78af_0 \\
        python3 amr_cassette_diagram.py --amrfinder-tsv ... --outfile ...

    # or with Singularity/Apptainer:
    singularity exec docker://quay.io/biocontainers/pygenomeviz:0.4.4--pyhdfd78af_0 \\
        python3 amr_cassette_diagram.py --amrfinder-tsv ... --outfile ...

Bandage panel note: `Bandage image` picks its output format from the file
extension. The training's own pipeline exports .svg (see run_workflow.sh),
which this script can't read directly -- either re-export as PNG:
    Bandage image assembly_graph.gfa assembly_graph.png
or convert the existing .svg, e.g. with cairosvg:
    pip install cairosvg
    cairosvg assembly_graph.svg -o assembly_graph.png

Usage:
    python3 amr_cassette_diagram.py \\
        --amrfinder-tsv ERR8282742.amrfinder.tsv \\
        --contig contig_1 \\
        --contig-length 6193746 \\
        --assembly-length 6808264 \\
        --bandage-image ERR8282742.flye_assembly_graph.png \\
        --outfile ERR8282742_amr_cassette.png
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.image import imread
from matplotlib.patches import Patch

from pygenomeviz import GenomeViz


# ---------------------------------------------------------------------
# Sizing -- tune these to make text/tracks bigger or smaller everywhere
# ---------------------------------------------------------------------
FIG_TRACK_HEIGHT = 2.2
FEATURE_TRACK_RATIO = 0.55
TRACK_LABELSIZE = 20
SUBLABEL_SIZE = 14
FEATURE_LABEL_SIZE = 13
SCALEBAR_LABELSIZE = 12
BACKBONE_LW = 4.0
FLANK_PAD_BP = 250            # bp of flanking sequence shown either side of the cassette
LEGEND_FONTSIZE = 12
CAPTION_FONTSIZE = 10
BANDAGE_CAPTION_FONTSIZE = 11
COMPOSITE_DPI = 200
BANDAGE_PANEL_HEIGHT = 4.0    # inches

# Fixed categorical colour order (validated for CVD-safe separation --
# BioCommons dataviz reference palette, slots 1-8). Classes are assigned
# slots in first-seen order along the cassette, so colours stay stable
# for a given TSV but aren't pinned to a specific class name -- edit
# CLASS_COLOUR_OVERRIDES below for a fixed, sample-independent mapping.
SLOT_COLOURS = [
    "#2a78d6",  # 1 blue
    "#eb6834",  # 2 orange
    "#1baf7a",  # 3 aqua
    "#eda100",  # 4 yellow
    "#e87ba4",  # 5 magenta
    "#008300",  # 6 green
    "#4a3aa7",  # 7 violet
    "#e34948",  # 8 red
]
CLASS_COLOUR_OVERRIDES: dict[str, str] = {}  # e.g. {"BETA-LACTAM": "#2a78d6"}

# AMRFinderPlus Method values meaning "alignment covers < 90% of the
# reference" (see the Method column docs) -- flagged in the gene label
# rather than a separate style channel, since opacity already carries
# %coverage continuously.
PARTIAL_METHODS = {"PARTIALX", "PARTIALP", "PARTIAL_CONTIG_ENDX", "PARTIAL_CONTIG_ENDP"}

MIN_ALPHA = 0.35
MAX_ALPHA = 1.0


@dataclass
class Hit:
    symbol: str
    name: str
    start: int
    stop: int
    strand: int
    element_class: str
    subclass: str
    method: str
    pct_coverage: float
    pct_identity: float
    ref_accession: str


def load_amrfinder_hits(tsv_path: str, contig: str | None) -> tuple[str, list[Hit]]:
    """Read an AMRFinderPlus TSV and return (contig_id_used, hits on it).

    If `contig` is None, auto-picks the contig with the most hits.
    """
    rows_by_contig: dict[str, list[dict]] = {}
    with open(tsv_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows_by_contig.setdefault(row["Contig id"], []).append(row)

    if not rows_by_contig:
        raise ValueError(f"No data rows found in {tsv_path}")

    if contig is None:
        contig = max(rows_by_contig, key=lambda c: len(rows_by_contig[c]))

    if contig not in rows_by_contig:
        raise ValueError(
            f"Contig '{contig}' not found in {tsv_path}. "
            f"Contigs present: {sorted(rows_by_contig)}"
        )

    hits = [
        Hit(
            symbol=row["Element symbol"],
            name=row["Element name"],
            start=int(row["Start"]),
            stop=int(row["Stop"]),
            strand=1 if row["Strand"] == "+" else -1,
            element_class=row["Class"],
            subclass=row["Subclass"],
            method=row["Method"],
            pct_coverage=float(row["% Coverage of reference"]),
            pct_identity=float(row["% Identity to reference"]),
            ref_accession=row["Closest reference accession"],
        )
        for row in rows_by_contig[contig]
    ]
    hits.sort(key=lambda h: h.start)
    return contig, hits


def assign_class_colours(hits: list[Hit]) -> dict[str, str]:
    """First-seen order -> fixed colour slots (unless overridden), so
    adjacent classes along the cassette get maximally distinct colours."""
    colours: dict[str, str] = {}
    for h in hits:
        if h.element_class in colours:
            continue
        colours[h.element_class] = CLASS_COLOUR_OVERRIDES.get(
            h.element_class, SLOT_COLOURS[len(colours) % len(SLOT_COLOURS)]
        )
    return colours


def coverage_to_alpha(pct_coverage: float) -> float:
    frac = max(0.0, min(100.0, pct_coverage)) / 100.0
    return MIN_ALPHA + frac * (MAX_ALPHA - MIN_ALPHA)


def render_cassette_panel(contig: str, hits: list[Hit], class_colours: dict[str, str],
                           outpath: Path) -> Path:
    lo = min(h.start for h in hits) - FLANK_PAD_BP
    hi = max(h.stop for h in hits) + FLANK_PAD_BP
    span = hi - lo
    cassette_span = max(h.stop for h in hits) - min(h.start for h in hits)

    # tick_style="bar" is pyGenomeViz 0.4.4's way of requesting a scale
    # bar (there is no separate set_scale_bar() method in this version).
    gv = GenomeViz(
        fig_track_height=FIG_TRACK_HEIGHT, feature_track_ratio=FEATURE_TRACK_RATIO,
        tick_style="bar", tick_labelsize=SCALEBAR_LABELSIZE,
    )
    track = gv.add_feature_track(
        contig, span, labelsize=TRACK_LABELSIZE,
        linewidth=BACKBONE_LW, linecolor="grey",
    )
    track.set_sublabel(
        f"{lo + 1:,}-{hi:,} bp excerpt (cassette span {cassette_span:,} bp)",
        size=SUBLABEL_SIZE,
    )

    for h in hits:
        label = h.symbol + (" (partial)" if h.method in PARTIAL_METHODS else "")
        track.add_feature(
            h.start - lo, h.stop - lo, h.strand,
            plotstyle="bigarrow",
            facecolor=class_colours[h.element_class], edgecolor="black", linewidth=0.6,
            label=label, labelsize=FEATURE_LABEL_SIZE, labelrotation=30,
            labelvpos="top", labelhpos="left",
            # alpha isn't a top-level add_feature() param in 0.4.4 -- pass
            # it through to the underlying matplotlib Patch instead.
            patch_kws={"alpha": coverage_to_alpha(h.pct_coverage)},
        )

    fig = gv.plotfig()
    fig.savefig(outpath, dpi=COMPOSITE_DPI, bbox_inches="tight")
    plt.close(fig)
    return outpath


def render_legend_strip(class_colours: dict[str, str], outpath: Path) -> Path:
    """Class-colour legend + the opacity caption, as their own small
    self-contained panel -- rendered separately (rather than layered onto
    the GenomeViz figure) so its spacing never collides with GenomeViz's
    own sublabel text."""
    fig, ax = plt.subplots(figsize=(9, 0.9))
    ax.axis("off")

    handles = [
        Patch(facecolor=c, edgecolor="black", label=cls.title())
        for cls, c in class_colours.items()
    ]
    ax.legend(
        handles=handles, loc="upper center", ncol=len(handles),
        frameon=False, fontsize=LEGEND_FONTSIZE, bbox_to_anchor=(0.5, 1.0),
    )
    ax.text(
        0.5, 0.05,
        "Arrow opacity = % coverage of the reference gene (faint = partial hit)",
        transform=ax.transAxes, ha="center", va="bottom",
        fontsize=CAPTION_FONTSIZE, color="dimgrey",
    )

    fig.savefig(outpath, dpi=COMPOSITE_DPI, bbox_inches="tight")
    plt.close(fig)
    return outpath


def compose_panels(image_heights: list[tuple[Path, float]], outfile: str) -> Path:
    """Stack pre-rendered panel images vertically (each already a clean,
    self-contained figure) into one final PNG."""
    n = len(image_heights)
    fig, axes = plt.subplots(
        n, 1, figsize=(9, sum(h for _, h in image_heights)),
        gridspec_kw={"height_ratios": [h for _, h in image_heights]},
    )
    axes = [axes] if n == 1 else list(axes)
    for ax, (img_path, _) in zip(axes, image_heights):
        ax.imshow(imread(img_path))
        ax.axis("off")

    fig.tight_layout()
    fig.savefig(outfile, dpi=COMPOSITE_DPI, bbox_inches="tight")
    plt.close(fig)
    return Path(outfile)


def render_bandage_panel(bandage_image: str, contig: str, contig_length: int | None,
                          assembly_length: int | None, outpath: Path) -> Path:
    fig, ax = plt.subplots(figsize=(9, BANDAGE_PANEL_HEIGHT))
    ax.imshow(imread(bandage_image))
    ax.axis("off")
    caption = f"Assembly graph (Bandage) -- {contig}"
    if contig_length and assembly_length:
        pct = 100 * contig_length / assembly_length
        caption += (
            f" is {contig_length:,} bp, {pct:.0f}% of the {assembly_length:,} bp "
            "assembly -- almost certainly the chromosome, not a plasmid."
        )
    ax.set_title(caption, fontsize=BANDAGE_CAPTION_FONTSIZE, loc="left")
    fig.savefig(outpath, dpi=COMPOSITE_DPI, bbox_inches="tight")
    plt.close(fig)
    return outpath


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--amrfinder-tsv", required=True, help="AMRFinderPlus TSV report")
    p.add_argument("--contig", default=None,
                    help="Contig id to draw (default: contig with the most hits)")
    p.add_argument("--bandage-image", default=None,
                    help="Pre-rendered Bandage PNG/JPG (see docstring for how to "
                         "get one) to show as a second panel")
    p.add_argument("--contig-length", type=int, default=None,
                    help="True length of the contig, for the Bandage-panel caption")
    p.add_argument("--assembly-length", type=int, default=None,
                    help="Total assembly length, for the Bandage-panel caption")
    p.add_argument("--outfile", default="amr_cassette.png")
    args = p.parse_args()

    contig, hits = load_amrfinder_hits(args.amrfinder_tsv, args.contig)
    if not hits:
        raise SystemExit(f"No AMRFinderPlus hits found on contig '{contig}'")

    class_colours = assign_class_colours(hits)

    outfile = Path(args.outfile)
    tmp_dir = outfile.parent
    cassette_png = tmp_dir / f"{outfile.stem}.cassette_panel.png"
    legend_png = tmp_dir / f"{outfile.stem}.legend_panel.png"
    render_cassette_panel(contig, hits, class_colours, cassette_png)
    render_legend_strip(class_colours, legend_png)

    panels = [(cassette_png, FIG_TRACK_HEIGHT), (legend_png, 0.9)]
    tmp_files = [cassette_png, legend_png]

    if args.bandage_image:
        bandage_png = tmp_dir / f"{outfile.stem}.bandage_panel.png"
        render_bandage_panel(
            args.bandage_image, contig, args.contig_length, args.assembly_length, bandage_png,
        )
        panels.append((bandage_png, BANDAGE_PANEL_HEIGHT))
        tmp_files.append(bandage_png)

    out = compose_panels(panels, str(outfile))
    for f in tmp_files:
        f.unlink(missing_ok=True)

    print(f"Saved {out}")
    print(f"Contig: {contig} | {len(hits)} hits | classes: {', '.join(class_colours)}")


if __name__ == "__main__":
    main()
