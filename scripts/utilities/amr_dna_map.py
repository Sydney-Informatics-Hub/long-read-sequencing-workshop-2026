#!/usr/bin/env python3
"""
amr_dna_map.py

Companion to amr_cassette_diagram.py: instead of base-pair-scale zoomed
cassette panels, draws a simple schematic "DNA map" for every contig in
an AMRFinderPlus TSV -- a thick ring for circular DNA, a straight
backbone line for non-circular DNA -- with each AMR/STRESS/VIRULENCE hit
drawn as a coloured segment at its true fractional position along the
molecule. Styling (drug-class colour, %coverage -> opacity, carbapenemase
red outline, core/plus scope linestyle) intentionally mirrors
amr_cassette_diagram.py's conventions, so the two scripts' output reads
as one consistent system -- but this script is otherwise fully
self-contained (no pyGenomeViz dependency: it isn't needed for a mock
circle/line, and skipping it means this script doesn't inherit that
script's pinned-to-0.4.4 constraint, see amr_cassette_diagram.py's
docstring). If you ever change a colour or style rule in one script,
change it in the other too.

The AMRFinderPlus TSV's "Contig id" is whatever the assembly FASTA it was
run on actually calls each sequence. If that assembly went through
merge_contigs.sh, the contigs were renamed (chromosome, plasmid_1,
plasmid_2, ...) -- names that appear nowhere in Flye's assembly_info.txt
(contig_1, contig_6, ...) or Plassembler's summary (its own numeric
plasmid ids), and Plassembler's numbering doesn't even line up with the
final plasmid_N suffixes, since merge_contigs.sh's numbering runs across
*all* plasmid-classified contigs (Flye's leftover circular repeats
included), not just Plassembler's. Pass that renaming with --id-map
(merge_contigs.sh writes it next to the draft assembly, as
"<draft_assembly>.id_map.txt") and every TSV contig id is translated back
to its raw Flye/Plassembler id before either file is consulted. Without
--id-map, TSV contig ids are assumed to already match the raw ids (true
only if the assembly was never renamed).

DNA length and circularity are then looked up by (translated) contig id,
checking two possible sources in turn:
    1. plassembler's *_summary.tsv (the final recovered plasmids, and the
       'chromosome' row plassembler carries through from the input flye
       assembly)
    2. flye's assembly_info.txt (the assembler's own raw contigs -- the
       only source for contigs plassembler didn't touch, e.g. repeat
       fragments an AMR hit happens to land on)
If a contig referenced in the AMRFinderPlus TSV isn't found in either
file, its length is inferred from the AMR hits' own coordinates (printed
as a warning) and it is assumed non-circular, since neither is knowable
without a real source.

One contig is treated as "the chromosome": either the one whose
AMRFinderPlus TSV contig id is literally 'chromosome' (as merge_contigs.sh
and plassembler's summary both name it), or -- if no contig has that id --
whichever resolved DNA molecule is both the largest one referenced in the
TSV and at least --chromosome-min-length bp (override either behaviour
with --chromosome-name). Every other contig
referenced in the TSV is drawn as a plasmid. The chromosome (if any hits
land on it) is written to its own PNG; every plasmid contig (if any) goes
into a second PNG as a small-multiples grid, one subplot per contig.

Usage:
    python3 amr_dna_map.py \\
        --amrfinder-tsv ERR8282752.amrfinder_150.tsv \\
        --flye-info flye/assembly_info.txt \\
        --plassembler-summary plassembler/plassembler_summary.tsv \\
        --id-map draft_assembly.fasta.id_map.txt \\
        --outdir .
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import imread, imsave
from matplotlib.patches import Wedge, Rectangle, Patch


# ---------------------------------------------------------------------
# Shared styling conventions, copied from amr_cassette_diagram.py (see
# the module docstring for why this is a copy rather than an import) --
# keep these in sync if either script's palette/encoding changes.
# ---------------------------------------------------------------------
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

CLASS_COLOUR_OVERRIDES: dict[str, str] = {
    "EFFLUX":         "#2a78d6",
    "BETA-LACTAM":    "#1baf7a",
    "AMINOGLYCOSIDE": "#eda100",
    "SULFONAMIDE":    "#008300",
    "QUINOLONE":      "#e87ba4",
    "FOSFOMYCIN":     "#8c5a2b",
    "PHENICOL":       "#2ba6c9",
    "METAL":     "#4a3aa7",
    "ACID":      "#eb6834",
    "VIRULENCE": "#e34948",
    "BIOCIDE":   "#7a3aa1",
    "HEAT":      "#8a9a1f",
}

PARTIAL_METHODS = {"PARTIALX", "PARTIALP", "PARTIAL_CONTIG_ENDX", "PARTIAL_CONTIG_ENDP"}
MIN_ALPHA = 0.35
MAX_ALPHA = 1.0

CARBAPENEMASE_SUBCLASS_KEYWORD = "CARBAPENEM"
CARBAPENEMASE_EDGECOLOR = "red"
CARBAPENEMASE_LINEWIDTH = 1.6
DEFAULT_EDGECOLOR = "black"
DEFAULT_LINEWIDTH = 1.6

CORE_SCOPE_LINESTYLE = "solid"
PLUS_SCOPE_LINESTYLE = (0, (6, 3))

# ---------------------------------------------------------------------
# This script's own sizing
# ---------------------------------------------------------------------
RING_WIDTH_FRAC = 0.18        # ring thickness as a fraction of the circle's radius
MIN_GENE_ANGLE_FRAC = 0.012   # minimum angular width (as a fraction of full circle) so a
                               # tiny/partial gene still shows up as a visible sliver
MIN_GENE_LINEAR_FRAC = 0.012  # same idea, for the linear backbone
MIN_GENE_LINEAR_GAP_FRAC = 0.006  # minimum gap between adjacent linear gene boxes,
                                    # so a tight cluster (e.g. an operon) doesn't
                                    # render as one overlapping, indistinguishable blob
LABEL_RADIUS_NEAR = 1.14      # circular label tiers (as a multiple of ring radius R):
LABEL_RADIUS_FAR = 1.32       # alternating near/far keeps adjacent labels from
                               # overlapping when hits sit close together (e.g. an operon)
MIN_GENE_LABEL_GAP_FRAC = 0.025  # minimum angular gap between adjacent circular labels,
                                   # as a fraction of the full circle (~9 degrees)
LABEL_FONTSIZE = 9
TITLE_FONTSIZE = 12
LEADER_LW = 0.8
BACKBONE_COLOUR = "#c9c9c9"
CIRCLE_FIGSIZE = (6.5, 6.5)
LINEAR_FIGSIZE = (10, 3)
PLASMID_SUBPLOT_SIZE = 5.0    # inches, per grid cell
PLASMID_GRID_MAX_COLS = 3
COMPOSITE_DPI = 200
CROP_PAD_PX = 15
GAP_IN = 0.15
LEGEND_FONTSIZE = 12
CAPTION_FONTSIZE = 10


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
    scope: str


@dataclass
class DNA:
    name: str
    length: int
    circular: bool
    source: str  # "plassembler", "flye", or "inferred"


# ---------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------
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
                scope=row["Scope"],
            )
            for row in rows
        ]
        hits.sort(key=lambda h: h.start)
        hits_by_contig[contig] = hits
    return hits_by_contig


def load_flye_info(path: str) -> dict[str, DNA]:
    out: dict[str, DNA] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            name = row["#seq_name"]
            out[name] = DNA(
                name=name,
                length=int(row["length"]),
                circular=row["circ."].strip().upper() == "Y",
                source="flye",
            )
    return out


def load_plassembler_summary(path: str) -> dict[str, DNA]:
    out: dict[str, DNA] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            name = row["contig"]
            out[name] = DNA(
                name=name,
                length=int(float(row["length"])),
                circular=row["circularity"].strip().lower() == "circular",
                source="plassembler",
            )
    return out


def load_id_map(path: str) -> dict[str, str]:
    """Read merge_contigs.sh's <draft_assembly>.id_map.txt (old_id<TAB>new_id,
    no header) and return the reverse mapping: new_id (the renamed contig id
    that actually appears in the AMRFinderPlus TSV, e.g. 'chromosome',
    'plasmid_1') -> old_id (the raw Flye/Plassembler contig id used in
    --flye-info / --plassembler-summary)."""
    out: dict[str, str] = {}
    with open(path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row:
                continue
            old_id, new_id = row[0], row[1]
            out[new_id] = old_id
    return out


def resolve_dna(contig_id: str, hits: list[Hit],
                 plassembler: dict[str, DNA],
                 flye: dict[str, DNA],
                 id_map: dict[str, str]) -> DNA:
    """plassembler's summary wins when a contig id appears in both (it's
    the more specific, post-reassembly source); flye's assembly_info is
    the fallback for anything plassembler didn't touch. contig_id is first
    translated to its raw Flye/Plassembler id via --id-map, since it may be
    a merge_contigs.sh-renamed id (chromosome, plasmid_1, ...) that appears
    in neither source file directly."""
    raw_id = id_map.get(contig_id, contig_id)
    if raw_id in plassembler:
        return plassembler[raw_id]
    if raw_id in flye:
        return flye[raw_id]
    inferred_length = max(h.stop for h in hits)
    raw_note = f" (raw id '{raw_id}')" if raw_id != contig_id else ""
    print(
        f"Warning: contig '{contig_id}'{raw_note} not found in --plassembler-summary "
        f"or --flye-info -- inferring length {inferred_length:,} bp from this "
        "contig's own AMR hit coordinates and assuming it is non-circular.",
        file=sys.stderr,
    )
    return DNA(name=contig_id, length=inferred_length, circular=False, source="inferred")


def classify_dna(contigs: dict[str, DNA], chromosome_name: str | None,
                  chromosome_min_length: int) -> tuple[str | None, list[str]]:
    if chromosome_name is not None:
        if chromosome_name not in contigs:
            raise SystemExit(
                f"--chromosome-name '{chromosome_name}' has no AMR/STRESS/VIRULENCE hits "
                f"in this TSV. Contigs with hits: {sorted(contigs)}"
            )
        chrom = chromosome_name
    else:
        named = [n for n in contigs if n.lower() == "chromosome"]
        if named:
            chrom = named[0]
        else:
            candidates = [n for n, r in contigs.items() if r.length >= chromosome_min_length]
            chrom = max(candidates, key=lambda n: contigs[n].length) if candidates else None
    plasmids = [n for n in contigs if n != chrom]
    return chrom, sorted(plasmids, key=lambda n: -contigs[n].length)


# ---------------------------------------------------------------------
# Per-hit styling (mirrors amr_cassette_diagram.py)
# ---------------------------------------------------------------------
def colour_group_key(h: Hit) -> str:
    return h.element_class if h.elem_type == "AMR" else h.subtype


def assign_class_colours(hits: list[Hit]) -> dict[str, str]:
    colours: dict[str, str] = {}
    for h in hits:
        key = colour_group_key(h)
        if key in colours:
            continue
        colours[key] = CLASS_COLOUR_OVERRIDES.get(
            key, SLOT_COLOURS[len(colours) % len(SLOT_COLOURS)]
        )
    return colours


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
    carbapenemase = is_carbapenemase(h)
    return {
        "edgecolor": CARBAPENEMASE_EDGECOLOR if carbapenemase else DEFAULT_EDGECOLOR,
        "linewidth": CARBAPENEMASE_LINEWIDTH if carbapenemase else DEFAULT_LINEWIDTH,
        "alpha": coverage_to_alpha(h.pct_coverage),
        "linestyle": CORE_SCOPE_LINESTYLE if h.scope == "core" else PLUS_SCOPE_LINESTYLE,
    }


# ---------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------
def draw_circular_dna(ax, dna: DNA, hits: list[Hit],
                       class_colours: dict[str, str]) -> None:
    """A thick ring standing in for the circular molecule (the 'wider
    edge' the two strands would occupy at this scale), with each hit
    recolouring the true angular slice of the ring it actually occupies --
    0 bp at 12 o'clock, running clockwise."""
    R = 1.0
    width = RING_WIDTH_FRAC * R
    ax.add_patch(Wedge((0, 0), R, 0, 360, width=width,
                        facecolor=BACKBONE_COLOUR, edgecolor="none", zorder=1))

    # Label angular positions are spread apart independently of the wedges
    # they belong to: a tight cluster of hits (e.g. an operon) has near-
    # identical true angles, so alternating radius tiers alone isn't enough
    # to keep their labels legible -- this nudges each label's angle just
    # far enough from its neighbours, while the wedge itself stays at the
    # true position and the leader line bends to connect the two.
    true_mid_fracs = [
        ((h.start / dna.length) + (h.stop / dna.length)) / 2 for h in hits
    ]
    label_mid_fracs = [
        x0 for x0, _ in _spread_intervals(
            [(m, m) for m in true_mid_fracs], MIN_GENE_LABEL_GAP_FRAC,
        )
    ]

    far = False
    for h, mid_frac, label_mid_frac in zip(hits, true_mid_fracs, label_mid_fracs):
        start_frac = h.start / dna.length
        stop_frac = h.stop / dna.length
        span_frac = max(stop_frac - start_frac, MIN_GENE_ANGLE_FRAC)
        angle_a = 90 - 360 * (mid_frac - span_frac / 2)
        angle_b = 90 - 360 * (mid_frac + span_frac / 2)
        theta1, theta2 = sorted((angle_a, angle_b))

        style = feature_style_kwargs(h)
        ax.add_patch(Wedge(
            (0, 0), R, theta1, theta2, width=width,
            facecolor=class_colours[colour_group_key(h)],
            edgecolor=style["edgecolor"], linewidth=style["linewidth"],
            alpha=style["alpha"], linestyle=style["linestyle"], zorder=2,
        ))

        anchor_angle_rad = math.radians(90 - 360 * mid_frac)
        anchor_x, anchor_y = R * math.cos(anchor_angle_rad), R * math.sin(anchor_angle_rad)
        label_angle_deg = 90 - 360 * label_mid_frac
        label_angle_rad = math.radians(label_angle_deg)
        label_r = LABEL_RADIUS_FAR if far else LABEL_RADIUS_NEAR
        far = not far
        label_x, label_y = label_r * math.cos(label_angle_rad), label_r * math.sin(label_angle_rad)
        ax.plot([anchor_x, label_x], [anchor_y, label_y], color="grey", lw=LEADER_LW, zorder=1)

        normalised = label_angle_deg % 360
        ha = "left" if normalised <= 90 or normalised >= 270 else "right"
        rot = label_angle_deg if ha == "left" else label_angle_deg + 180
        ax.text(label_x, label_y, gene_label(h), fontsize=LABEL_FONTSIZE,
                rotation=rot, rotation_mode="anchor", ha=ha, va="center", zorder=3)

    ax.set_xlim(-2.1, 2.1)
    ax.set_ylim(-2.1, 2.1)
    ax.set_aspect("equal")
    ax.axis("off")


def _spread_intervals(spans: list[tuple[float, float]], min_gap: float) -> list[tuple[float, float]]:
    """Nudge a start-sorted list of (x0, x1) spans rightward just enough that
    none overlap or sit closer than min_gap apart. Without this, a tight
    cluster of hits (e.g. an operon) whose minimum-width boxes overlap draws
    as a single indistinguishable blob -- this keeps each one individually
    legible while preserving order and true relative position as closely as
    the minimum spacing allows."""
    out: list[tuple[float, float]] = []
    prev_x1 = -math.inf
    for x0, x1 in spans:
        if x0 < prev_x1 + min_gap:
            shift = (prev_x1 + min_gap) - x0
            x0, x1 = x0 + shift, x1 + shift
        out.append((x0, x1))
        prev_x1 = x1
    return out


def draw_linear_dna(ax, dna: DNA, hits: list[Hit],
                     class_colours: dict[str, str]) -> None:
    """A straight backbone standing in for a non-circular contig, with
    each hit drawn as a coloured block at its true position along it."""
    ax.plot([0, dna.length], [0, 0], color=BACKBONE_COLOUR,
            linewidth=14, solid_capstyle="butt", zorder=1)

    min_span = dna.length * MIN_GENE_LINEAR_FRAC
    min_gap = dna.length * MIN_GENE_LINEAR_GAP_FRAC
    raw_spans = []
    for h in hits:
        span = max(h.stop - h.start, min_span)
        mid = (h.start + h.stop) / 2
        raw_spans.append((mid - span / 2, mid + span / 2))
    spans = _spread_intervals(raw_spans, min_gap)

    above = True
    for h, (x0, x1) in zip(hits, spans):
        style = feature_style_kwargs(h)
        ax.add_patch(Rectangle(
            (x0, -7), x1 - x0, 14,
            facecolor=class_colours[colour_group_key(h)],
            edgecolor=style["edgecolor"], linewidth=style["linewidth"],
            alpha=style["alpha"], linestyle=style["linestyle"], zorder=2,
        ))

        mid = (x0 + x1) / 2
        label_y = 16 if above else -16
        va = "bottom" if above else "top"
        ax.plot([mid, mid], [7 if above else -7, label_y], color="grey", lw=LEADER_LW, zorder=1)
        ax.text(mid, label_y, gene_label(h), fontsize=LABEL_FONTSIZE,
                rotation=90, ha="center", va=va, zorder=3)
        above = not above

    pad = max(dna.length * 0.05, 1)
    rightmost = max(dna.length, spans[-1][1] if spans else dna.length)
    ax.set_xlim(-pad, rightmost + pad)
    ax.set_ylim(-45, 45)
    ax.axis("off")


def draw_dna(ax, name: str, dna: DNA, hits: list[Hit],
             class_colours: dict[str, str], label_suffix: str = "") -> None:
    if dna.circular:
        draw_circular_dna(ax, dna, hits, class_colours)
    else:
        draw_linear_dna(ax, dna, hits, class_colours)
    topology = "circular" if dna.circular else "linear"
    subtitle = (
        f"{topology}, {dna.length:,} bp, {len(hits)} hit{'s' if len(hits) != 1 else ''}"
    )
    if dna.source == "inferred":
        subtitle += " (length inferred from hits -- not in either source file)"
    ax.set_title(f"{name}{label_suffix}\n{subtitle}", fontsize=TITLE_FONTSIZE)


# ---------------------------------------------------------------------
# Legend + compositing (mirrors amr_cassette_diagram.py's helpers)
# ---------------------------------------------------------------------
def render_legend_strip(class_colours: dict[str, str], hits: list[Hit], outpath: Path) -> Path:
    fig, ax = plt.subplots(figsize=(CIRCLE_FIGSIZE[0], 2.35))
    ax.axis("off")

    type_by_key: dict[str, str] = {}
    for h in hits:
        type_by_key.setdefault(colour_group_key(h), h.elem_type)

    def _class_patch(cls: str, c: str) -> Patch:
        return Patch(facecolor=c, edgecolor="black",
                      label=f"{type_by_key.get(cls, '').upper()}: "
                            f"{'Other/Unclassified' if cls == 'NA' else cls.title()}")

    # class_colours is the global palette (shared across the chromosome and
    # plasmid figures so colours match between them), but `hits`/type_by_key
    # is scoped to just this one figure -- so both branches must also check
    # a key is actually present in type_by_key, not just what type it is,
    # or a class from the *other* figure's contigs leaks into this legend.
    amr_handles = [_class_patch(cls, c) for cls, c in class_colours.items()
                   if type_by_key.get(cls) == "AMR"]
    other_handles = [_class_patch(cls, c) for cls, c in class_colours.items()
                      if cls in type_by_key and type_by_key.get(cls) != "AMR"]

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

    if any(is_carbapenemase(h) for h in hits) and "BETA-LACTAM" in class_colours:
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
        "Ring/backbone position is a true fraction of DNA length, not to a fixed bp scale",
    ]
    for i, caption in enumerate(captions):
        ax.text(
            0.5, 0.2 - i * 0.09, caption,
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
    imgs = [_autocrop(_to_rgba_uint8(imread(p))) for p in image_paths]
    canvas_w = max(im.shape[1] for im in imgs)
    gap_px = round(GAP_IN * COMPOSITE_DPI)
    gap_row = np.full((gap_px, canvas_w, 4), 255, dtype=np.uint8)

    rows = []
    for i, im in enumerate(imgs):
        h, w = im.shape[:2]
        if w < canvas_w:
            left_pad = (canvas_w - w) // 2
            right_pad = canvas_w - w - left_pad
            im = np.concatenate([
                np.full((h, left_pad, 4), 255, dtype=np.uint8), im,
                np.full((h, right_pad, 4), 255, dtype=np.uint8),
            ], axis=1)
        if i > 0:
            rows.append(gap_row)
        rows.append(im)

    imsave(outfile, np.concatenate(rows, axis=0))
    return Path(outfile)


# ---------------------------------------------------------------------
# Top-level rendering
# ---------------------------------------------------------------------
def render_chromosome(name: str, dna: DNA, hits: list[Hit],
                       class_colours: dict[str, str], outpath: Path) -> Path:
    figsize = CIRCLE_FIGSIZE if dna.circular else LINEAR_FIGSIZE
    fig, ax = plt.subplots(figsize=figsize)
    draw_dna(ax, name, dna, hits, class_colours, label_suffix=" (chromosome)")

    tmp_dir = outpath.parent
    main_png = tmp_dir / f"{outpath.stem}.main.png"
    legend_png = tmp_dir / f"{outpath.stem}.legend.png"
    fig.savefig(main_png, dpi=COMPOSITE_DPI, bbox_inches="tight")
    plt.close(fig)
    render_legend_strip(class_colours, hits, legend_png)

    out = compose_panels([main_png, legend_png], str(outpath))
    main_png.unlink(missing_ok=True)
    legend_png.unlink(missing_ok=True)
    return out


def render_plasmids(plasmid_names: list[str], contigs: dict[str, DNA],
                     hits_by_contig: dict[str, list[Hit]], class_colours: dict[str, str],
                     outpath: Path) -> Path:
    n = len(plasmid_names)
    ncols = min(PLASMID_GRID_MAX_COLS, n)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(PLASMID_SUBPLOT_SIZE * ncols, PLASMID_SUBPLOT_SIZE * nrows),
        squeeze=False,
    )
    for i, name in enumerate(plasmid_names):
        ax = axes[i // ncols][i % ncols]
        draw_dna(ax, name, contigs[name], hits_by_contig[name], class_colours)
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.tight_layout()

    tmp_dir = outpath.parent
    main_png = tmp_dir / f"{outpath.stem}.main.png"
    legend_png = tmp_dir / f"{outpath.stem}.legend.png"
    fig.savefig(main_png, dpi=COMPOSITE_DPI, bbox_inches="tight")
    plt.close(fig)

    plasmid_hits = [h for name in plasmid_names for h in hits_by_contig[name]]
    render_legend_strip(class_colours, plasmid_hits, legend_png)

    out = compose_panels([main_png, legend_png], str(outpath))
    main_png.unlink(missing_ok=True)
    legend_png.unlink(missing_ok=True)
    return out


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--amrfinder-tsv", required=True, help="AMRFinderPlus TSV report")
    p.add_argument("--flye-info", default=None,
                    help="flye's assembly_info.txt (for contig length/circularity)")
    p.add_argument("--plassembler-summary", default=None,
                    help="plassembler's *_summary.tsv (for contig length/circularity; "
                         "takes precedence over --flye-info when a contig id is in both)")
    p.add_argument("--id-map", default=None,
                    help="merge_contigs.sh's <draft_assembly>.id_map.txt (old_id<TAB>new_id, "
                         "no header). Required whenever the AMRFinderPlus TSV's contig ids "
                         "(e.g. chromosome, plasmid_1) were renamed by merge_contigs.sh and so "
                         "differ from the raw ids in --flye-info / --plassembler-summary.")
    p.add_argument("--chromosome-name", default=None,
                    help="Force this contig id to be treated as the chromosome (default: "
                         "a contig literally named 'chromosome', else the largest contig "
                         "referenced in the TSV that is over --chromosome-min-length)")
    p.add_argument("--chromosome-min-length", type=int, default=1_000_000,
                    help="Minimum length in bp for a contig to be auto-classified as the "
                         "chromosome when none is named 'chromosome' (default: 1,000,000)")
    p.add_argument("--outdir", default=".",
                    help="Directory to write the chromosome/plasmid map PNGs into")
    args = p.parse_args()

    if not args.flye_info and not args.plassembler_summary:
        raise SystemExit(
            "Provide at least one of --flye-info / --plassembler-summary, otherwise every "
            "contig's length has to be inferred from its own AMR hits, which is unreliable."
        )

    hits_by_contig = load_all_contigs(args.amrfinder_tsv)
    plassembler = load_plassembler_summary(args.plassembler_summary) if args.plassembler_summary else {}
    flye = load_flye_info(args.flye_info) if args.flye_info else {}
    id_map = load_id_map(args.id_map) if args.id_map else {}

    contigs: dict[str, DNA] = {
        contig_id: resolve_dna(contig_id, hits, plassembler, flye, id_map)
        for contig_id, hits in hits_by_contig.items()
    }

    chrom_name, plasmid_names = classify_dna(
        contigs, args.chromosome_name, args.chromosome_min_length,
    )

    all_hits = [h for hits in hits_by_contig.values() for h in hits]
    class_colours = assign_class_colours(all_hits)

    tsv_path = Path(args.amrfinder_tsv)
    tsv_stem = tsv_path.name[:-len(".tsv")] if tsv_path.name.endswith(".tsv") else tsv_path.stem
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if chrom_name is not None:
        chrom_png = outdir / f"{tsv_stem}.chromosome_map.png"
        out = render_chromosome(
            chrom_name, contigs[chrom_name], hits_by_contig[chrom_name], class_colours, chrom_png,
        )
        print(f"Saved {out}  ({chrom_name}: {len(hits_by_contig[chrom_name])} hits)")
    else:
        print("No chromosome-sized contig found among the AMR-hit contigs -- skipping the "
              "chromosome map. Use --chromosome-name to force one, or lower "
              "--chromosome-min-length.")

    if plasmid_names:
        plasmid_png = outdir / f"{tsv_stem}.plasmid_map.png"
        out = render_plasmids(plasmid_names, contigs, hits_by_contig, class_colours, plasmid_png)
        print(f"Saved {out}  ({len(plasmid_names)} contig(s): {', '.join(plasmid_names)})")
    else:
        print("No plasmid-classified contigs found among the AMR-hit contigs -- skipping "
              "the plasmid map.")


if __name__ == "__main__":
    main()
