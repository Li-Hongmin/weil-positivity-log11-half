#!/bin/sh
# Typeset paper/paper.md into paper/paper.pdf with pandoc + tectonic.
#
# paper.md is the manuscript exactly as submitted and is not modified.
# The title and abstract are taken from metadata.json (the submission
# metadata; the manuscript body carries no abstract section). The body
# is read from the first "## " heading onward, because the title and the
# author line are supplied through the YAML header. Two wide displays,
# (6.2) and (7.1), are typeset at a smaller size so that they fit the
# text width; their content is unchanged.
set -eu
cd "$(dirname "$0")"
mkdir -p build
python3 - <<'PY'
import json, re
m = json.load(open("metadata.json"))
def q(s): return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'
front = [
    "---",
    "title: " + q(m["title"]),
    'author: "Hongmin Li (Institute of Science Tokyo; The University of Tokyo)"',
    "abstract: " + q(m["abstract"]),
    'date: "rStation version 2, November 2026"',
    "documentclass: article",
    "fontsize: 11pt",
    "geometry: margin=1in",
    "header-includes: |",
    "  \\usepackage{amsmath,amssymb}",
    "  \\usepackage{xurl}",
    "  \\allowdisplaybreaks",
    "---",
    "",
]
text = open("paper.md", encoding="utf-8").read()
body = text[text.index("\n## ") + 1:]

def shrink(match):
    inner = match.group(1)
    if "\\begin{array}{c|cc}" in inner:
        size = "tiny"
    elif "73.817062678870728" in inner:
        size = "scriptsize"
    else:
        return match.group(0)
    return "\\begin{%s}\n\\[%s\\]\n\\end{%s}" % (size, inner, size)

body = re.sub(r"^\$\$\n(.*?)\n\$\$$", shrink, body, flags=re.S | re.M)
open("build/paper-typeset.md", "w", encoding="utf-8").write("\n".join(front) + body)
PY
pandoc build/paper-typeset.md -f markdown -t pdf --pdf-engine=tectonic -o paper.pdf
rm -rf build
echo "wrote paper.pdf"
