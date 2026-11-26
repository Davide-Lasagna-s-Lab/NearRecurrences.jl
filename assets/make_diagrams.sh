#!/usr/bin/env bash
# Draws the diagrams of the README, assets/<name>.svg, from their TikZ sources in
# assets/diagrams/<name>.tex, by tectonic and pdftocairo (poppler).
#
#     bash assets/make_diagrams.sh
set -euo pipefail
cd "$(dirname "$0")/diagrams"
for tex in *.tex; do
    name="${tex%.tex}"
    tectonic --chatter minimal "$tex"
    pdftocairo -svg "$name.pdf" "../$name.svg"
    rm "$name.pdf"
done
