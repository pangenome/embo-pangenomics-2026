#!/usr/bin/env bash
set -euo pipefail

mkdir -p graphs figures results

if command -v odgi >/dev/null 2>&1; then
    ODGI_BIN="$(command -v odgi)"
else
    ODGI_BIN="$(find /gnu/store -path '*/bin/odgi' -type f 2>/dev/null | head -1 || true)"
fi

if [[ -z "${ODGI_BIN}" ]]; then
    echo "Could not find odgi. Try: guix install odgi" >&2
    exit 1
fi

echo "Using odgi: ${ODGI_BIN}" >&2
printf "graph\tlength\tnodes\tedges\tpaths\tsteps\tpng\n" > results/odgi_graph_stats.tsv

for gfa in "$@"; do
    [[ -s "$gfa" ]] || continue
    base="$(basename "$gfa" .gfa)"
    og="graphs/${base}.og"
    png="figures/${base}.odgi.png"

    echo "odgi build/viz ${gfa}" >&2
    "$ODGI_BIN" build -g "$gfa" -o "$og" 2> "results/${base}.odgi_build.log"
    "$ODGI_BIN" viz -i "$og" -o "$png" 2> "results/${base}.odgi_viz.log"
    "$ODGI_BIN" stats -i "$og" -S \
    | awk -v graph="$gfa" -v png="$png" 'NR == 2 {print graph "\t" $0 "\t" png}' \
    >> results/odgi_graph_stats.tsv
done
