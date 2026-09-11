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
    # Default: every contig in the TSV gets its own PNG, auto-named next
    # to the input TSV (<tsv-name-without-.tsv>.<contig>.png):
    python3 amr_cassette_diagram.py --amrfinder-tsv ERR8282742.amrfinder.tsv

    # One specific contig, with a Bandage panel and a chosen output name
    # (--outfile only works together with --contig -- it names one PNG):
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
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import imread, imsave
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
FLANK_PAD_BP = 250            # bp of flanking sequence shown either side of a zoomed cluster
LEGEND_FONTSIZE = 12
CAPTION_FONTSIZE = 10
BANDAGE_CAPTION_FONTSIZE = 11
COMPOSITE_DPI = 200
BANDAGE_PANEL_HEIGHT = 4.0    # inches
CANVAS_WIDTH_IN = 9.0         # shared figure width -- literally every panel (overview,
                               # legend, bandage, every zoomed cluster) renders at exactly
                               # this width. Earlier this script instead shrank narrow
                               # zoomed panels to their proportional width, which kept the
                               # bp-per-inch scale true but ran into a pyGenomeViz 0.4.4
                               # rendering bug: a bigarrow's arrowhead is sized from a fixed
                               # fraction of *panel width in inches*, not from the gene's own
                               # length, so a narrow panel gets a near-zero-length,
                               # full-height arrowhead -- a degenerate sliver that reads as a
                               # stray spike poking above/below the backbone (confirmed by
                               # rendering the same gene at several widths: the spike shrinks
                               # to negligible once the panel is back at full canvas width).
                               # Fix: keep every panel at this same full width and instead
                               # widen the *view* around small loci (see display_span in
                               # render_cassette_panel) so bp-per-inch stays identical
                               # everywhere without ever shrinking a panel.
SINGLE_GENE_TRACK_HEIGHT = 1.3  # shorter row for a lone-gene panel (only one label to fit)
SCALE_REF_TRACK_HEIGHT = 1.0    # standalone scale-bar panel -- see render_scale_reference_panel
GAP_IN = 0.15                  # blank vertical gap between stacked panels in the composite
CROP_PAD_PX = 15               # breathing room left after autocropping a panel -- see _autocrop

# Real AMR "cassettes" (integrons, resistance islands, co-located operons)
# are physically clustered genes. AMRFinderPlus TSVs, though, often also
# carry unrelated single-gene hits scattered across the rest of the
# replicon (a chromosomal efflux pump here, a metal-resistance regulator
# there). Plotting everything on one track sized to min(start)..max(stop)
# squeezes the real cluster into an illegibly small sliver of pixels to
# make room for those far-apart singletons. So hits are grouped by gap
# first: every location -- cluster or lone gene -- gets its own zoomed
# panel at true scale, and a whole-contig overview ties them together.
CLUSTER_GAP_BP = 15_000        # hits within this many bp of each other count as one cluster
MIN_CLUSTER_SIZE_FOR_ZOOM = 1  # every cluster gets a zoomed panel, lone genes included
OVERVIEW_COLLAPSE_MIN_SIZE = 2  # multi-gene clusters collapse to a box on the overview;
                                 # lone genes still show as their own small coloured arrow
OVERVIEW_FLANK_BP = 20_000     # bp of padding either side when no --contig-length is given
OVERVIEW_TRACK_HEIGHT = 1.6
OVERVIEW_FEATURE_LABEL_SIZE = 11
OVERVIEW_MIN_LABEL_GAP_PX = 45  # markers rendered closer than this in the overview get
                                 # merged into one combined marker -- otherwise their
                                 # (now-vertical) labels have nowhere to go but on top of
                                 # each other. Zoom panels are unaffected; this only changes
                                 # how the overview groups markers for display.

# Fixed categorical colour order (validated for CVD-safe separation --
# BioCommons dataviz reference palette, slots 1-8), used as a fallback for
# any colour_group_key CLASS_COLOUR_OVERRIDES doesn't cover -- assigned in
# first-seen order along the cassette, so an uncovered key still gets a
# consistent colour *within* one run, just not guaranteed across runs.
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

# Fixed, sample-independent colour per colour_group_key (an AMR Class, or
# a STRESS/VIRULENCE Subtype) -- every key seen across the ERR8282742/751/
# 752/753 reports is pinned here, so e.g. AMR/BETA-LACTAM is always the
# same aqua whether it's the only hit on a contig or one of twenty, and
# whichever sample it's in. Add a new key here (rather than letting it
# fall through to SLOT_COLOURS) the first time a new Class/Subtype shows
# up, to keep that guarantee as more samples get run.
CLASS_COLOUR_OVERRIDES: dict[str, str] = {
    # AMR, coloured by Class
    "EFFLUX":         "#2a78d6",  # blue
    "BETA-LACTAM":    "#1baf7a",  # aqua
    "AMINOGLYCOSIDE": "#eda100",  # gold
    "SULFONAMIDE":    "#008300",  # green
    "QUINOLONE":      "#e87ba4",  # magenta
    "FOSFOMYCIN":     "#8c5a2b",  # brown
    "PHENICOL":       "#2ba6c9",  # sky

    # STRESS / VIRULENCE, coloured by Subtype
    "METAL":     "#4a3aa7",  # violet
    "ACID":      "#eb6834",  # orange
    "VIRULENCE": "#e34948",  # red
    "BIOCIDE":   "#7a3aa1",  # purple
    "HEAT":      "#8a9a1f",  # olive
}

# AMRFinderPlus Method values meaning "alignment covers < 90% of the
# reference" (see the Method column docs) -- flagged in the gene label
# rather than a separate style channel, since opacity already carries
# %coverage continuously.
PARTIAL_METHODS = {"PARTIALX", "PARTIALP", "PARTIAL_CONTIG_ENDX", "PARTIAL_CONTIG_ENDP"}

MIN_ALPHA = 0.35
MAX_ALPHA = 1.0

# Carbapenemases are the single highest-priority AMR finding clinically
# (they defeat last-line carbapenem therapy), so they get their own
# outline colour independent of drug-class fill colour: a gene reads as
# "carbapenemase" at a glance regardless of which class colour it happens
# to have. Detected via Subclass rather than a fixed gene-symbol list, so
# any carbapenemase family (IMP, NDM, KPC, OXA-48-like, VIM, ...) is caught.
CARBAPENEMASE_SUBCLASS_KEYWORD = "CARBAPENEM"
CARBAPENEMASE_EDGECOLOR = "red"
CARBAPENEMASE_LINEWIDTH = 3.2
DEFAULT_EDGECOLOR = "black"
DEFAULT_LINEWIDTH = 1.6

# AMRFinderPlus "Scope" column: "core" is the curated, high-confidence AMR
# gene set; "plus" is the broader/less-curated set (stress, biocide, extra
# virulence genes). Encoded as outline style -- solid vs a bold, widely
# spaced dash -- rather than colour, since colour already carries drug
# class and red is already spoken for by carbapenemase status. (A dashed
# outline on a narrow bigarrow used to spike past the tip -- see the
# CANVAS_WIDTH_IN note above for why every panel is now full-width, which
# fixes this too.)
CORE_SCOPE_LINESTYLE = "solid"
PLUS_SCOPE_LINESTYLE = (0, (6, 3))


@dataclass
class Hit:
    symbol: str
    name: str
    start: int
    stop: int
    strand: int
    elem_type: str
    subtype: str
    element_class: str
    subclass: str
    method: str
    pct_coverage: float
    pct_identity: float
    ref_accession: str
    scope: str


def load_all_contigs(tsv_path: str) -> dict[str, list[Hit]]:
    """Read an AMRFinderPlus TSV and return every contig's hits, keyed by
    contig id, each sorted by start position."""
    rows_by_contig: dict[str, list[dict]] = {}
    with open(tsv_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows_by_contig.setdefault(row["Contig id"], []).append(row)

    if not rows_by_contig:
        raise ValueError(f"No data rows found in {tsv_path}")

    hits_by_contig: dict[str, list[Hit]] = {}
    for contig, rows in rows_by_contig.items():
        hits = [
            Hit(
                symbol=row["Element symbol"],
                name=row["Element name"],
                start=int(row["Start"]),
                stop=int(row["Stop"]),
                strand=1 if row["Strand"] == "+" else -1,
                elem_type=row["Type"],
                subtype=row["Subtype"],
                element_class=row["Class"],
                subclass=row["Subclass"],
                method=row["Method"],
                pct_coverage=float(row["% Coverage of reference"]),
                pct_identity=float(row["% Identity to reference"]),
                ref_accession=row["Closest reference accession"],
                scope=row["Scope"],
            )
            for row in rows
        ]
        hits.sort(key=lambda h: h.start)
        hits_by_contig[contig] = hits
    return hits_by_contig


def _natural_sort_key(s: str) -> list:
    """'contig_2' before 'contig_10' -- plain string sort gets that
    backwards, since '1' < '10' < '2' lexicographically."""
    return [int(tok) if tok.isdigit() else tok.lower() for tok in re.split(r"(\d+)", s)]


def colour_group_key(h: Hit) -> str:
    """What a gene's colour is keyed on. AMR-typed hits get one colour per
    drug Class (AMINOGLYCOSIDE, BETA-LACTAM, EFFLUX, ...) -- that's the
    informative axis for an AMR gene. Everything else (STRESS, VIRULENCE)
    gets one colour per Subtype instead of per Class: STRESS/METAL hits in
    particular don't need a separate colour for every individual metal
    (ARSENIC vs COPPER vs MERCURY, ...) -- one "Metal" colour reads better,
    and the gene symbol/label already says which metal. This also means
    the same Class curated under two Types doesn't collide: EFFLUX genes
    AMRFinderPlus files under AMR (mexE, acrF, ...) get the AMR/Class
    colour, while the STRESS/BIOCIDE EFFLUX genes (ttgA/B/R) get the
    STRESS/BIOCIDE colour instead, since they're grouped by Subtype."""
    return h.element_class if h.elem_type == "AMR" else h.subtype


def assign_class_colours(hits: list[Hit]) -> dict[str, str]:
    """First-seen order -> fixed colour slots (unless overridden), so
    adjacent classes along the cassette get maximally distinct colours."""
    colours: dict[str, str] = {}
    for h in hits:
        key = colour_group_key(h)
        if key in colours:
            continue
        colours[key] = CLASS_COLOUR_OVERRIDES.get(
            key, SLOT_COLOURS[len(colours) % len(SLOT_COLOURS)]
        )
    return colours


def cluster_hits(hits: list[Hit], max_gap: int) -> list[list[Hit]]:
    """Group hits (already sorted by start) into clusters, starting a new
    cluster whenever a hit starts more than `max_gap` bp past the furthest
    stop seen so far in the current cluster."""
    clusters: list[list[Hit]] = []
    cluster_end = None
    for h in hits:
        if clusters and h.start - cluster_end <= max_gap:
            clusters[-1].append(h)
            cluster_end = max(cluster_end, h.stop)
        else:
            clusters.append([h])
            cluster_end = h.stop
    return clusters


def group_for_overview(clusters: list[list[Hit]], span: int) -> list[list[list[Hit]]]:
    """Merge clusters that would render closer together than
    OVERVIEW_MIN_LABEL_GAP_PX in the overview into one combined group, so
    their (vertical) labels don't collide -- this only changes how the
    overview groups markers for display; zoom panels are built from the
    original `clusters`, unaffected."""
    min_gap_bp = span * OVERVIEW_MIN_LABEL_GAP_PX / (CANVAS_WIDTH_IN * COMPOSITE_DPI)
    groups: list[list[list[Hit]]] = []
    group_end = None
    for c in clusters:
        c_start = min(h.start for h in c)
        c_stop = max(h.stop for h in c)
        if groups and c_start - group_end <= min_gap_bp:
            groups[-1].append(c)
            group_end = max(group_end, c_stop)
        else:
            groups.append([c])
            group_end = c_stop
    return groups


def _format_index_range(idxs: list[int]) -> str:
    if len(idxs) > 1 and idxs == list(range(idxs[0], idxs[0] + len(idxs))):
        return f"{idxs[0]}-{idxs[-1]}"
    return ",".join(str(i) for i in idxs)


def coverage_to_alpha(pct_coverage: float) -> float:
    frac = max(0.0, min(100.0, pct_coverage)) / 100.0
    return MIN_ALPHA + frac * (MAX_ALPHA - MIN_ALPHA)


def is_carbapenemase(h: Hit) -> bool:
    return CARBAPENEMASE_SUBCLASS_KEYWORD in h.subclass.upper()


def gene_label(h: Hit) -> str:
    label = h.symbol
    if h.method in PARTIAL_METHODS:
        label += " (partial)"
    return label


def feature_style_kwargs(h: Hit) -> dict:
    """edgecolor/linewidth/patch_kws shared by every add_feature() call for
    a Hit, so the carbapenemase and core/plus encodings stay consistent
    between the zoomed panels and the overview's singleton arrows."""
    carbapenemase = is_carbapenemase(h)
    linestyle = CORE_SCOPE_LINESTYLE if h.scope == "core" else PLUS_SCOPE_LINESTYLE
    return {
        "edgecolor": CARBAPENEMASE_EDGECOLOR if carbapenemase else DEFAULT_EDGECOLOR,
        "linewidth": CARBAPENEMASE_LINEWIDTH if carbapenemase else DEFAULT_LINEWIDTH,
        "patch_kws": {"alpha": coverage_to_alpha(h.pct_coverage), "linestyle": linestyle},
    }


def render_cassette_panel(contig: str, hits: list[Hit], class_colours: dict[str, str],
                           outpath: Path, track_label: str | None = None,
                           show_scale_bar: bool = True,
                           display_span: int | None = None,
                           fig_track_height: float = FIG_TRACK_HEIGHT) -> Path:
    gene_lo = min(h.start for h in hits)
    gene_hi = max(h.stop for h in hits)
    cassette_span = gene_hi - gene_lo

    if display_span is None:
        # No shared scale to match (e.g. the single-cluster case) --
        # just pad by a fixed flank either side of the genes themselves.
        lo, hi = gene_lo - FLANK_PAD_BP, gene_hi + FLANK_PAD_BP
    else:
        # Every zoomed panel is centered in a window of the *same* size
        # (see CANVAS_WIDTH_IN) so bp-per-inch is identical everywhere,
        # whether this cluster fills the window or sits as a small blob
        # in the middle of a lot of flanking sequence.
        center = (gene_lo + gene_hi) / 2
        lo = round(center - display_span / 2)
        hi = lo + display_span
    span = hi - lo

    # tick_style="bar" is pyGenomeViz 0.4.4's way of requesting a scale
    # bar (there is no separate set_scale_bar() method in this version).
    # Only one panel (the widest) needs show_scale_bar=True: since every
    # panel shares the same scale, repeating the bar on each is clutter.
    gv = GenomeViz(
        fig_width=CANVAS_WIDTH_IN, fig_track_height=fig_track_height,
        feature_track_ratio=FEATURE_TRACK_RATIO,
        tick_style="bar" if show_scale_bar else None, tick_labelsize=SCALEBAR_LABELSIZE,
    )
    track = gv.add_feature_track(
        track_label or contig, span, labelsize=TRACK_LABELSIZE,
        linewidth=BACKBONE_LW, linecolor="grey",
    )
    track.set_sublabel(
        f"{lo + 1:,}-{hi:,} bp excerpt (cluster span {cassette_span:,} bp)",
        size=SUBLABEL_SIZE,
    )

    for h in hits:
        track.add_feature(
            h.start - lo, h.stop - lo, h.strand,
            plotstyle="bigarrow",
            facecolor=class_colours[colour_group_key(h)],
            label=gene_label(h), labelsize=FEATURE_LABEL_SIZE, labelrotation=90,
            labelvpos="top", labelhpos="left",
            **feature_style_kwargs(h),
        )

    fig = gv.plotfig()
    fig.savefig(outpath, dpi=COMPOSITE_DPI, bbox_inches="tight")
    plt.close(fig)
    return outpath


def render_overview_panel(contig: str, hits: list[Hit], clusters: list[list[Hit]],
                           zoom_index_by_id: dict[int, int], class_colours: dict[str, str],
                           contig_length: int | None, outpath: Path) -> Path:
    """Whole-contig context track. Sparse singleton hits are drawn as their
    own small arrows; any cluster dense enough to get its own zoomed panel
    is instead collapsed to one labelled box here, so the two panels
    together show both where things are and what they actually look like."""
    if contig_length:
        lo, hi = 0, contig_length
    else:
        lo = min(h.start for h in hits) - OVERVIEW_FLANK_BP
        hi = max(h.stop for h in hits) + OVERVIEW_FLANK_BP
    span = hi - lo

    gv = GenomeViz(
        fig_width=CANVAS_WIDTH_IN, fig_track_height=OVERVIEW_TRACK_HEIGHT,
        feature_track_ratio=FEATURE_TRACK_RATIO,
        tick_style="bar", tick_labelsize=SCALEBAR_LABELSIZE,
    )
    track = gv.add_feature_track(
        contig, span, labelsize=TRACK_LABELSIZE,
        linewidth=BACKBONE_LW, linecolor="grey",
    )
    track.set_sublabel(
        f"Whole-contig context -- {len(hits)} hits in {len(clusters)} location"
        f"{'s' if len(clusters) != 1 else ''} (all zoomed below)",
        size=SUBLABEL_SIZE,
    )

    # Two clusters can be genomically distinct (kept as separate zoom
    # panels below) but still sit only a handful of pixels apart at
    # whole-contig scale -- group_for_overview merges those into one
    # combined marker here so their vertical labels don't overlap.
    for group in group_for_overview(clusters, span):
        if len(group) == 1 and len(group[0]) < OVERVIEW_COLLAPSE_MIN_SIZE:
            # plotstyle="box", not "bigarrow": at whole-contig scale a
            # single gene is sub-pixel wide, so an arrowhead has no shaft
            # left to attach to -- it degenerates into a spike through the
            # backbone (the same underlying pyGenomeViz sizing quirk noted
            # under CANVAS_WIDTH_IN above, just unavoidable here since we
            # can't widen the whole-contig view the way the zoomed panels
            # do). Strand is meaningless at this scale anyway, so a plain
            # box -- like the collapsed multi-gene cluster markers below --
            # is both cleaner and more honest about what's visible here.
            cluster = group[0]
            h = cluster[0]
            zoom_idx = zoom_index_by_id[id(cluster)]
            label = f"[{zoom_idx}] " + gene_label(h)
            track.add_feature(
                h.start - lo, h.stop - lo, 1,
                plotstyle="box",
                facecolor=class_colours[colour_group_key(h)],
                label=label, labelsize=OVERVIEW_FEATURE_LABEL_SIZE, labelrotation=90,
                labelvpos="top", labelhpos="left",
                **feature_style_kwargs(h),
            )
        else:
            group_hits = [h for cluster in group for h in cluster]
            c_start = min(h.start for h in group_hits)
            c_stop = max(h.stop for h in group_hits)
            has_cp = any(is_carbapenemase(h) for h in group_hits)
            idxs = sorted(zoom_index_by_id[id(c)] for c in group)
            label = f"[{_format_index_range(idxs)}] {len(group_hits)} genes"
            if len(group) > 1:
                label += f" ({len(group)} locations)"
            if has_cp:
                label += " (incl. carbapenemase)"
            track.add_feature(
                c_start - lo, c_stop - lo, 1,
                plotstyle="box",
                facecolor="dimgrey",
                edgecolor=CARBAPENEMASE_EDGECOLOR if has_cp else DEFAULT_EDGECOLOR,
                linewidth=CARBAPENEMASE_LINEWIDTH if has_cp else 0.8,
                label=label,
                labelsize=OVERVIEW_FEATURE_LABEL_SIZE, labelrotation=90,
                labelvpos="top", labelhpos="left",
            )

    fig = gv.plotfig()
    fig.savefig(outpath, dpi=COMPOSITE_DPI, bbox_inches="tight")
    plt.close(fig)
    return outpath


def render_scale_reference_panel(display_span: int, outpath: Path) -> Path:
    """A standalone scale bar, not attached to any particular gene
    cluster's row -- putting the one shared scale bar on an arbitrary
    zoomed panel (the widest one, say) reads as if it belongs to that
    cluster specifically, not to every panel below it. This panel has no
    features at all, just a bare backbone at the shared span/width, so
    there's nothing to attribute the bar to but the scale itself."""
    gv = GenomeViz(
        fig_width=CANVAS_WIDTH_IN, fig_track_height=SCALE_REF_TRACK_HEIGHT,
        feature_track_ratio=FEATURE_TRACK_RATIO,
        tick_style="bar", tick_labelsize=SCALEBAR_LABELSIZE,
    )
    track = gv.add_feature_track(
        "Scale", display_span, labelsize=TRACK_LABELSIZE,
        linewidth=BACKBONE_LW, linecolor="grey",
    )
    track.set_sublabel(
        "Scale bar below applies to every zoomed panel that follows "
        "(all share one bp-per-inch scale)",
        size=SUBLABEL_SIZE,
    )
    fig = gv.plotfig()
    fig.savefig(outpath, dpi=COMPOSITE_DPI, bbox_inches="tight")
    plt.close(fig)
    return outpath


def render_legend_strip(class_colours: dict[str, str], hits: list[Hit], outpath: Path) -> Path:
    """Class-colour legend + the encoding captions, as their own small
    self-contained panel -- rendered separately (rather than layered onto
    the GenomeViz figure) so its spacing never collides with GenomeViz's
    own sublabel text."""
    fig, ax = plt.subplots(figsize=(CANVAS_WIDTH_IN, 2.35))
    ax.axis("off")

    # Each colour key is either a Class (AMR hits) or a Subtype (everything
    # else) -- see colour_group_key. Label it with its Type too (AMR /
    # STRESS / VIRULENCE) so the legend shows which axis each swatch came
    # from, not just the bare class/subtype name.
    type_by_key: dict[str, str] = {}
    for h in hits:
        type_by_key.setdefault(colour_group_key(h), h.elem_type)

    def _class_patch(cls: str, c: str) -> Patch:
        return Patch(facecolor=c, edgecolor="black",
                      label=f"{type_by_key.get(cls, '').upper()}: "
                            f"{'Other/Unclassified' if cls == 'NA' else cls.title()}")

    # AMR gets its own row, everything else (STRESS, VIRULENCE, ...) gets a
    # second row below it -- two clearly separate axes (what drug beats
    # this gene vs. what non-drug pressure it survives) read better apart
    # than interleaved along one line. Each row is however wide it needs
    # to be: with enough classes (or long labels like "Copper/Silver")
    # that can exceed CANVAS_WIDTH_IN. That's fine -- an Axes clips its
    # legend to its own bounding box by default, which silently truncated
    # swatches here before set_clip_on(False); with clipping off, each row
    # renders at its full natural width and bbox_inches='tight' captures
    # all of it, so _autocrop and compose_panels' pad-to-widest-panel
    # logic can do their job on the real content instead of a fragment.
    amr_handles = [_class_patch(cls, c) for cls, c in class_colours.items()
                   if type_by_key.get(cls) == "AMR"]
    other_handles = [_class_patch(cls, c) for cls, c in class_colours.items()
                      if type_by_key.get(cls) != "AMR"]

    row_y = 1.0
    for handles in (amr_handles, other_handles):
        if not handles:
            continue
        row_legend = ax.legend(
            handles=handles, loc="upper center", ncol=len(handles),
            frameon=False, fontsize=LEGEND_FONTSIZE, bbox_to_anchor=(0.5, row_y),
        )
        row_legend.set_clip_on(False)
        ax.add_artist(row_legend)
        row_y -= 0.14

    # Carbapenem isn't a class of its own -- it's a *subclass* of
    # BETA-LACTAM, so it gets its own row below the class rows (rather
    # than sitting inline in one of them) filled with the Beta-Lactam
    # colour, red outline on top, to show it's a flagged subset of that
    # class.
    if any(is_carbapenemase(h) for h in hits):
        carbapenem_legend = ax.legend(
            handles=[Patch(
                facecolor=class_colours["BETA-LACTAM"], edgecolor=CARBAPENEMASE_EDGECOLOR,
                linewidth=CARBAPENEMASE_LINEWIDTH, label="Carbapenem (subset)",
            )],
            loc="upper center", ncol=1,
            frameon=False, fontsize=LEGEND_FONTSIZE, bbox_to_anchor=(0.5, 0.6),
        )
        carbapenem_legend.set_clip_on(False)
        ax.add_artist(carbapenem_legend)

    # Scope is a separate, orthogonal axis from class -- the same class
    # can contain both 'core' and 'plus' hits (e.g. BETA-LACTAM genes are
    # sometimes one, sometimes the other, depending on the specific gene
    # matched), so it can't be folded into the class rows above without
    # duplicating classes across a "core" and "plus" line. A row of
    # generic swatches shows the outline style directly instead of only
    # describing it in the caption text below.
    scope_handles = [
        Patch(facecolor="white", edgecolor=DEFAULT_EDGECOLOR, linewidth=DEFAULT_LINEWIDTH,
              linestyle=CORE_SCOPE_LINESTYLE, label="Core scope"),
        Patch(facecolor="white", edgecolor=DEFAULT_EDGECOLOR, linewidth=DEFAULT_LINEWIDTH,
              linestyle=PLUS_SCOPE_LINESTYLE, label="Plus scope"),
    ]
    scope_legend = ax.legend(
        handles=scope_handles, loc="upper center", ncol=2,
        frameon=False, fontsize=LEGEND_FONTSIZE, bbox_to_anchor=(0.5, 0.44),
    )
    scope_legend.set_clip_on(False)

    captions = [
        "Arrow opacity = % coverage of the reference gene (faint = partial hit)",
        "Red outline = carbapenem subclass (Class BETA-LACTAM, Subclass contains CARBAPENEM)",
        "All zoomed panels share one bp-per-inch scale, shown once as its own "
        "scale-bar panel above",
    ]
    for i, caption in enumerate(captions):
        ax.text(
            0.5, 0.2 - i * 0.09,
            caption,
            transform=ax.transAxes, ha="center", va="bottom",
            fontsize=CAPTION_FONTSIZE, color="dimgrey",
        )

    fig.savefig(outpath, dpi=COMPOSITE_DPI, bbox_inches="tight")
    plt.close(fig)
    return outpath


def _to_rgba_uint8(arr: np.ndarray) -> np.ndarray:
    if arr.dtype != np.uint8:
        arr = (arr * 255).round().astype(np.uint8)
    if arr.ndim == 2:
        arr = np.stack([arr] * 3, axis=-1)
    if arr.shape[-1] == 3:
        alpha = np.full((*arr.shape[:2], 1), 255, dtype=np.uint8)
        arr = np.concatenate([arr, alpha], axis=-1)
    return arr


def _autocrop(arr: np.ndarray, pad_px: int = CROP_PAD_PX) -> np.ndarray:
    """Trim near-white margin down to a small fixed pad on all four sides.
    bbox_inches='tight' is not actually tight for these panels -- pyGenomeViz
    reserves vertical space per track (e.g. for the tick/scale-bar area, or
    room a rotated label *could* need) well beyond what it draws, and a
    legend forced into one row can leave most of its declared width blank.
    That reserved-but-empty space is what made the composite so tall;
    cropping to real ink fixes it in one place instead of chasing each
    panel type's own padding quirk."""
    non_white = np.any(arr[:, :, :3] < 250, axis=-1)
    rows = np.where(non_white.any(axis=1))[0]
    cols = np.where(non_white.any(axis=0))[0]
    if rows.size == 0 or cols.size == 0:
        return arr
    h, w = arr.shape[:2]
    top, bottom = max(0, rows[0] - pad_px), min(h, rows[-1] + 1 + pad_px)
    left, right = max(0, cols[0] - pad_px), min(w, cols[-1] + 1 + pad_px)
    return arr[top:bottom, left:right]


def compose_panels(image_paths: list[Path], outfile: str) -> Path:
    """Stack pre-rendered panel images vertically into one final PNG,
    pixel-for-pixel -- left-aligned on a shared canvas width rather than
    each stretched to fill it, and autocropped first (see _autocrop)."""
    imgs = [_autocrop(_to_rgba_uint8(imread(p))) for p in image_paths]
    canvas_w = max(im.shape[1] for im in imgs)
    gap_px = round(GAP_IN * COMPOSITE_DPI)
    gap_row = np.full((gap_px, canvas_w, 4), 255, dtype=np.uint8)

    rows = []
    for i, im in enumerate(imgs):
        h, w = im.shape[:2]
        if w < canvas_w:
            pad = np.full((h, canvas_w - w, 4), 255, dtype=np.uint8)
            im = np.concatenate([im, pad], axis=1)
        if i > 0:
            rows.append(gap_row)
        rows.append(im)

    imsave(outfile, np.concatenate(rows, axis=0))
    return Path(outfile)


def render_bandage_panel(bandage_image: str, contig: str, contig_length: int | None,
                          assembly_length: int | None, outpath: Path) -> Path:
    fig, ax = plt.subplots(figsize=(CANVAS_WIDTH_IN, BANDAGE_PANEL_HEIGHT))
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


def render_amr_cassette(contig: str, hits: list[Hit], outfile: Path,
                         bandage_image: str | None, contig_length: int | None,
                         assembly_length: int | None) -> Path:
    """Build every panel for one contig and compose them into `outfile`."""
    class_colours = assign_class_colours(hits)
    clusters = cluster_hits(hits, CLUSTER_GAP_BP)
    zoom_clusters = [c for c in clusters if len(c) >= MIN_CLUSTER_SIZE_FOR_ZOOM]

    tmp_dir = outfile.parent
    legend_png = tmp_dir / f"{outfile.stem}.legend_panel.png"
    render_legend_strip(class_colours, hits, legend_png)
    panels = []
    tmp_files = [legend_png]

    def cluster_span(cluster: list[Hit]) -> int:
        return (max(h.stop for h in cluster) - min(h.start for h in cluster)) \
            + 2 * FLANK_PAD_BP

    if len(clusters) == 1:
        # Every hit is already one physical cluster -- no scattered
        # singletons to separate out, so a whole-contig overview would
        # just duplicate the zoomed panel.
        cassette_png = tmp_dir / f"{outfile.stem}.cassette_panel.png"
        render_cassette_panel(contig, hits, class_colours, cassette_png)
        panels.append(cassette_png)
        tmp_files.append(cassette_png)
    else:
        zoom_index_by_id = {id(c): i for i, c in enumerate(zoom_clusters, start=1)}

        overview_png = tmp_dir / f"{outfile.stem}.overview_panel.png"
        render_overview_panel(
            contig, hits, clusters, zoom_index_by_id, class_colours,
            contig_length, overview_png,
        )
        panels.append(overview_png)
        tmp_files.append(overview_png)

        # One shared bp-per-inch scale across every zoomed panel: every
        # panel is rendered at the same CANVAS_WIDTH_IN, each centered on
        # its own cluster but showing the same total bp span as the
        # widest cluster needs -- so a lone gene appears as a small blob
        # in a lot of flanking sequence rather than the panel itself
        # shrinking (see the CANVAS_WIDTH_IN comment for why: a shrunk
        # panel triggers a pyGenomeViz rendering bug). The scale bar is
        # shown once, on its own standalone panel -- putting it on one
        # arbitrary cluster's row (even the widest) reads as if it
        # belongs to that cluster specifically, not to every panel below.
        display_span = max(cluster_span(c) for c in zoom_clusters)

        scale_png = tmp_dir / f"{outfile.stem}.scale_panel.png"
        render_scale_reference_panel(display_span, scale_png)
        panels.append(scale_png)
        tmp_files.append(scale_png)

        for cluster in zoom_clusters:
            i = zoom_index_by_id[id(cluster)]
            zoom_png = tmp_dir / f"{outfile.stem}.cluster_{i}_panel.png"
            track_label = f"{contig} -- [{i}/{len(zoom_clusters)}]"
            track_height = SINGLE_GENE_TRACK_HEIGHT if len(cluster) == 1 else FIG_TRACK_HEIGHT
            render_cassette_panel(
                contig, cluster, class_colours, zoom_png, track_label,
                show_scale_bar=False,
                display_span=display_span, fig_track_height=track_height,
            )
            panels.append(zoom_png)
            tmp_files.append(zoom_png)

    panels.append(legend_png)

    if bandage_image:
        bandage_png = tmp_dir / f"{outfile.stem}.bandage_panel.png"
        render_bandage_panel(bandage_image, contig, contig_length, assembly_length, bandage_png)
        panels.append(bandage_png)
        tmp_files.append(bandage_png)

    out = compose_panels(panels, str(outfile))
    for f in tmp_files:
        f.unlink(missing_ok=True)

    print(f"Saved {out}")
    print(
        f"Contig: {contig} | {len(hits)} hits in {len(clusters)} location(s), "
        f"{len(zoom_clusters)} zoomed | classes: {', '.join(class_colours)}"
    )
    return out


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--amrfinder-tsv", required=True, help="AMRFinderPlus TSV report")
    p.add_argument("--contig", default=None,
                    help="Contig id to draw (default: every contig in the TSV, each its "
                         "own PNG)")
    p.add_argument("--bandage-image", default=None,
                    help="Pre-rendered Bandage PNG/JPG (see docstring for how to "
                         "get one) to show as a second panel")
    p.add_argument("--contig-length", type=int, default=None,
                    help="True length of the contig, for the Bandage-panel caption")
    p.add_argument("--assembly-length", type=int, default=None,
                    help="Total assembly length, for the Bandage-panel caption")
    p.add_argument("--outfile", default=None,
                    help="Output PNG path. Only valid together with --contig (a single "
                         "contig); with no --contig, every contig gets its own "
                         "auto-named PNG next to the input TSV (<tsv-name-without-.tsv>."
                         "<contig>.png) and --outfile must be omitted.")
    args = p.parse_args()

    hits_by_contig = load_all_contigs(args.amrfinder_tsv)

    if args.contig is not None:
        if args.contig not in hits_by_contig:
            raise SystemExit(
                f"Contig '{args.contig}' not found in {args.amrfinder_tsv}. "
                f"Contigs present: {sorted(hits_by_contig)}"
            )
        contigs = [args.contig]
    else:
        if args.outfile is not None:
            raise SystemExit(
                "--outfile requires --contig (it names a single PNG). Omit --outfile "
                "to render every contig, each to its own auto-named PNG, or pass "
                "--contig to render just one."
            )
        contigs = sorted(hits_by_contig, key=_natural_sort_key)

    tsv_path = Path(args.amrfinder_tsv)
    tsv_stem = tsv_path.name[:-len(".tsv")] if tsv_path.name.endswith(".tsv") else tsv_path.stem

    for contig in contigs:
        hits = hits_by_contig[contig]
        if not hits:
            continue
        outfile = Path(args.outfile) if args.outfile else tsv_path.parent / f"{tsv_stem}.{contig}.png"
        render_amr_cassette(
            contig, hits, outfile,
            args.bandage_image, args.contig_length, args.assembly_length,
        )


if __name__ == "__main__":
    main()
