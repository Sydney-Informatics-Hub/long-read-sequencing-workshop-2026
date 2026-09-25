#!/usr/bin/env python3
"""
filter_assembly_graph.py

Filter a Flye or Plassembler GFA (v1) assembly graph to a set of contig IDs.
Uses only the Python standard library (Python >= 3.7).

Flye
  - Graph:   assembly_graph.gfa   (segments are named edge_N)
  - Mapping: assembly_info.txt    (graph_path column maps contig_N / scaffold_N
                                   to signed edge IDs; '*' and '??' are gaps)
  - Contig -> edge resolution order:
      1. assembly_info.txt, if supplied with --info
      2. P (path) lines in the GFA, if present (Flye >= 2.9)
      3. naming convention contig_N -> edge_N
  - IDs given directly as edge_N are used as segment names.

Plassembler
  - Graph:   <prefix>_plasmids.gfa  (Unicycler graph; segment names are the
                                     contig numbers used in _plasmids.fasta)
  - Summary: <prefix>_summary.tsv   (optional; used only to check that the
                                     requested IDs are listed)

Links (L), containments (C) and paths (P) are kept only when every segment
they reference is retained. Header (H) lines are always kept.

Examples
  python filter_assembly_graph.py -g flye_out/assembly_graph.gfa \\
      -i flye_out/assembly_info.txt --ids contig_3,contig_7 -o sub.gfa

  python filter_assembly_graph.py -g plassembler_plasmids.gfa \\
      -s plassembler_summary.tsv --ids-file keep.txt -o plasmids_sub.gfa

  # also keep segments up to 1 link away from the selected segments
  python filter_assembly_graph.py -g assembly_graph.gfa --ids contig_2 \\
      --neighbours 1 -o sub.gfa
"""

import argparse
import csv
import sys
from collections import defaultdict, deque

GAP_TOKENS = {"*", "??"}


# ----------------------------------------------------------------- input ----

def read_ids(ids_arg, ids_file):
    ids = []
    if ids_arg:
        ids.extend(x.strip() for x in ids_arg.split(","))
    if ids_file:
        with open(ids_file) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                # accept whitespace/comma separated, or a FASTA-style header
                for tok in line.replace(",", " ").split():
                    ids.append(tok.lstrip(">"))
    ids = [i for i in ids if i]
    # preserve order, remove duplicates
    return list(dict.fromkeys(ids))


def read_gfa(path):
    """Return (lines, segment_names, adjacency, paths)."""
    lines = []
    segments = set()
    adjacency = defaultdict(set)
    paths = {}
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            lines.append(line)
            f = line.split("\t")
            rec = f[0]
            if rec == "S" and len(f) >= 2:
                segments.add(f[1])
            elif rec in ("L", "C") and len(f) >= 5:
                a, b = f[1], f[3]
                adjacency[a].add(b)
                adjacency[b].add(a)
            elif rec == "P" and len(f) >= 3:
                paths[f[1]] = [s[:-1] for s in f[2].split(",") if s[-1] in "+-"]
    return lines, segments, adjacency, paths


def read_flye_info(path):
    """Map contig/scaffold name -> list of edge_N segment names."""
    mapping = {}
    with open(path) as fh:
        header = None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not f or not f[0]:
                continue
            if f[0].startswith("#"):
                header = [h.lstrip("#").strip() for h in f]
                continue
            if header and "graph_path" in header:
                gp = f[header.index("graph_path")]
            else:
                gp = f[-1]  # graph_path is the last column in Flye 2.x
            edges = []
            for tok in gp.split(","):
                tok = tok.strip()
                if not tok or tok in GAP_TOKENS:
                    continue
                edges.append("edge_" + tok.lstrip("+-"))
            mapping[f[0]] = edges
    return mapping


def read_plassembler_summary(path):
    """Return the set of contig IDs listed in a Plassembler summary TSV."""
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header is None:
            return set()
        lower = [h.strip().lower() for h in header]
        col = lower.index("contig") if "contig" in lower else 0
        return {row[col].strip() for row in reader if len(row) > col}


# ------------------------------------------------------------ resolution ----

def detect_format(segments):
    if segments and all(s.startswith("edge_") for s in segments):
        return "flye"
    return "plassembler"


def resolve_flye(ids, segments, info_map, gfa_paths):
    keep, missing, report = set(), [], []
    for cid in ids:
        if cid in segments:
            edges, source = [cid], "segment name"
        elif info_map is not None and cid in info_map:
            edges, source = info_map[cid], "assembly_info.txt"
        elif cid in gfa_paths:
            edges, source = gfa_paths[cid], "GFA P line"
        elif cid.startswith("contig_") and "edge_" + cid[7:] in segments:
            edges, source = ["edge_" + cid[7:]], "contig_N -> edge_N convention"
        else:
            missing.append(cid)
            continue
        absent = [e for e in edges if e not in segments]
        present = list(dict.fromkeys(e for e in edges if e in segments))
        keep.update(present)
        report.append((cid, present, source, absent))
    return keep, missing, report


def resolve_plassembler(ids, segments, summary_ids):
    keep, missing, report = set(), [], []
    for cid in ids:
        if cid in segments:
            keep.add(cid)
            note = []
            if summary_ids is not None and cid not in summary_ids:
                note = ["not listed in summary TSV"]
            report.append((cid, [cid], "segment name", note))
        else:
            missing.append(cid)
    return keep, missing, report


def expand_neighbours(seeds, adjacency, depth):
    if depth <= 0:
        return set(seeds)
    seen = set(seeds)
    queue = deque((s, 0) for s in seeds)
    while queue:
        node, d = queue.popleft()
        if d == depth:
            continue
        for nb in adjacency.get(node, ()):
            if nb not in seen:
                seen.add(nb)
                queue.append((nb, d + 1))
    return seen


# ---------------------------------------------------------------- output ----

def filter_lines(lines, keep):
    out, dropped = [], defaultdict(int)
    for line in lines:
        f = line.split("\t")
        rec = f[0]
        if rec == "H":
            out.append(line)
        elif rec == "S":
            if f[1] in keep:
                out.append(line)
            else:
                dropped["S"] += 1
        elif rec in ("L", "C"):
            if f[1] in keep and f[3] in keep:
                out.append(line)
            else:
                dropped[rec] += 1
        elif rec == "P":
            segs = [s[:-1] for s in f[2].split(",") if s and s[-1] in "+-"]
            if segs and all(s in keep for s in segs):
                out.append(line)
            else:
                dropped["P"] += 1
        else:
            dropped[rec] += 1
    return out, dropped


def main():
    ap = argparse.ArgumentParser(
        description="Filter a Flye or Plassembler GFA graph to selected contigs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Examples", 1)[1] if "Examples" in __doc__ else None,
    )
    ap.add_argument("-g", "--graph", required=True, help="input GFA file")
    ap.add_argument("-o", "--output", required=True, help="output GFA file")
    ap.add_argument("-f", "--format", choices=["auto", "flye", "plassembler"],
                    default="auto", help="assembler that produced the graph (default: auto)")
    ap.add_argument("-i", "--info", help="Flye assembly_info.txt")
    ap.add_argument("-s", "--summary", help="Plassembler _summary.tsv")
    ids_grp = ap.add_argument_group("contig IDs (at least one required)")
    ids_grp.add_argument("--ids", help="comma-separated contig IDs")
    ids_grp.add_argument("--ids-file", help="file with one or more IDs per line")
    ap.add_argument("-n", "--neighbours", type=int, default=0,
                    help="also keep segments within N links of the selection (default: 0)")
    ap.add_argument("--strict", action="store_true",
                    help="exit with an error if any requested ID is not found")
    args = ap.parse_args()

    ids = read_ids(args.ids, args.ids_file)
    if not ids:
        ap.error("no contig IDs supplied (use --ids and/or --ids-file)")

    lines, segments, adjacency, gfa_paths = read_gfa(args.graph)
    if not segments:
        sys.exit("ERROR: no S (segment) lines found in " + args.graph)

    fmt = detect_format(segments) if args.format == "auto" else args.format
    log = sys.stderr
    print(f"Format: {fmt}  |  segments in input: {len(segments)}", file=log)

    if fmt == "flye":
        info_map = read_flye_info(args.info) if args.info else None
        keep, missing, report = resolve_flye(ids, segments, info_map, gfa_paths)
    else:
        summary_ids = read_plassembler_summary(args.summary) if args.summary else None
        keep, missing, report = resolve_plassembler(ids, segments, summary_ids)

    for cid, segs, source, notes in report:
        msg = f"  {cid} -> {','.join(segs) if segs else '(none)'}  [{source}]"
        if notes:
            label = "missing edges: " if fmt == "flye" else ""
            msg += f"  WARNING {label}{','.join(notes)}"
        print(msg, file=log)
    if missing:
        print(f"WARNING: IDs not found: {', '.join(missing)}", file=log)
        if args.strict:
            sys.exit(1)
    if not keep:
        sys.exit("ERROR: no segments selected; output not written.")

    selected = len(keep)
    keep = expand_neighbours(keep, adjacency, args.neighbours)
    if args.neighbours > 0:
        print(f"Neighbour expansion (depth {args.neighbours}): "
              f"{selected} -> {len(keep)} segments", file=log)

    out, dropped = filter_lines(lines, keep)
    with open(args.output, "w") as fh:
        fh.write("\n".join(out) + "\n")

    kept_counts = defaultdict(int)
    for line in out:
        kept_counts[line.split("\t", 1)[0]] += 1
    summary = ", ".join(f"{k}={v}" for k, v in sorted(kept_counts.items()))
    removed = ", ".join(f"{k}={v}" for k, v in sorted(dropped.items())) or "none"
    print(f"Written {args.output}  |  kept: {summary}  |  removed: {removed}", file=log)


if __name__ == "__main__":
    main()
