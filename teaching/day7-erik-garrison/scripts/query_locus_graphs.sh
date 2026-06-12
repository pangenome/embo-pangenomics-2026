#!/usr/bin/env bash
set -euo pipefail

mkdir -p graphs results

syng="${SYN_G:-yeast235.syng}"
agc="${AGC:-yeast235.agc}"

run_graph() {
    local label="$1"
    local region="$2"
    local engine="$3"
    local suffix="$4"
    shift 4

    local out="graphs/${label}.${suffix}.gfa"
    if [[ -s "$out" ]]; then
        echo "skipping existing ${out}" >&2
        return
    fi
    echo "building ${out}" >&2
    impg query \
        -a "$syng" \
        --sequence-files "$agc" \
        -r "$region" \
        -d 1k \
        -o "gfa:${engine}" \
        "$@" \
    > "$out"
}

# Controls: one core gene and three pangenome-flavored loci.
run_graph ADE2  "SGDref#0#chrXV:564476-566191"  poa  "poa.ext0"    --syng-extension 0
run_graph ADE2  "SGDref#0#chrXV:564476-566191"  poa  "poa.ext1000" --syng-extension 1000
run_graph ADE2  "SGDref#0#chrXV:564476-566191"  syng "syng.ext1000" --syng-extension 1000
run_graph ADE2  "SGDref#0#chrXV:564476-566191"  pggb "pggb.ext1000" --syng-extension 1000

run_graph SUC2  "SGDref#0#chrIX:36000-40000" poa  "poa.ext1000" --syng-extension 1000
run_graph SUC2  "SGDref#0#chrIX:36000-40000" syng "syng.ext1000" --syng-extension 1000
run_graph SUC2  "SGDref#0#chrIX:36000-40000" pggb "pggb.ext1000" --syng-extension 1000

run_graph CUP1  "SGDref#0#chrVIII:210000-216000" syng "syng.ext1000" --syng-extension 1000
run_graph CUP1  "SGDref#0#chrVIII:210000-216000" syng "syng.ext0" --syng-extension 0
run_graph CUP1  "SGDref#0#chrVIII:210000-216000" syng "syng.ext5000" --syng-extension 5000
run_graph CUP1  "SGDref#0#chrVIII:210000-216000" syng "syng.seed_drop0" --syng-extension 1000 --syng-seed-drop-top-fraction 0
run_graph CUP1  "SGDref#0#chrVIII:210000-216000" syng "syng.seed_drop01" --syng-extension 1000 --syng-seed-drop-top-fraction 0.01
run_graph CUP1  "SGDref#0#chrVIII:210000-216000" syng "syng.walk3" --syng-extension 1000 --syng-seed-walk-anchors 3
run_graph CUP1  "SGDref#0#chrVIII:210000-216000" syng "syng.walk7" --syng-extension 1000 --syng-seed-walk-anchors 7

run_graph FLO11 "SGDref#0#chrIX:386000-397000" syng "syng.ext1000" --syng-extension 1000

run_graph ENA1  "SGDref#0#chrIV:532000-540000" syng "syng.ext1000" --syng-extension 1000

run_graph HXT67 "SGDref#0#chrIV:1152000-1163000" syng "syng.ext1000" --syng-extension 1000

{
    printf "graph\tsegments\tlinks\tpaths\twalks\tbytes\n"
    for gfa in graphs/*.gfa; do
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$gfa" \
            "$(grep -c '^S' "$gfa" || true)" \
            "$(grep -c '^L' "$gfa" || true)" \
            "$(grep -c '^P' "$gfa" || true)" \
            "$(grep -c '^W' "$gfa" || true)" \
            "$(wc -c < "$gfa")"
    done
} > results/locus_graph_gfa_counts.tsv
