#!/usr/bin/env python3
"""
Concatenate multiple SVG files into a single HTML file,
with a heading (the file's basename) above each SVG.

Usage:
    python svg_to_html.py output.html file1.svg file2.svg file3.svg ...
    python svg_to_html.py output.html *.svg
"""

import sys
import os


def build_html(svg_paths):
    parts = [
        "<!DOCTYPE html>",
        "<html>",
        "<head>",
        "<meta charset='utf-8'><title>SVG Gallery</title>",
        "<style>",
        "  .svg-container { width: 75vw; }",
        "  .svg-container svg { width: 100%; height: auto; display: block; }",
        "</style>",
        "</head>",
        "<body>",
    ]

    for path in svg_paths:
        basename = os.path.basename(path)
        with open(path, "r", encoding="utf-8") as f:
            svg_content = f.read()

        parts.append(f"<h2>{basename}</h2>")
        parts.append("<div class='svg-container'>")
        parts.append(svg_content)
        parts.append("</div>")
        parts.append("<hr>")

    parts.append("</body>")
    parts.append("</html>")

    return "\n".join(parts)


def main():
    if len(sys.argv) < 3:
        print("Usage: python svg_to_html.py output.html file1.svg file2.svg ...")
        sys.exit(1)

    output_path = sys.argv[1]
    svg_paths = sys.argv[2:]

    html = build_html(svg_paths)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"Wrote {len(svg_paths)} SVG(s) to {output_path}")


if __name__ == "__main__":
    main()
