#!/usr/bin/env bash
set -euo pipefail

chrom="${1:-SGDref#0#chrV}"
window="${2:-10000}"
distance="${3:-1k}"
syng="${4:-yeast235.syng}"
agc="${5:-yeast235.agc}"
fai="${6:-yeast235.fa.gz.fai}"
min_hit_len="${7:-5000}"

length=$(awk -v chrom="$chrom" '$1 == chrom {print $2}' "$fai")
if [[ -z "${length}" ]]; then
    echo "Could not find ${chrom} in ${fai}" >&2
    exit 1
fi

ref_chrom="${chrom##*#}"
canonical_chrom_re='^chr(I|II|III|IV|V|VI|VII|VIII|IX|X|XI|XII|XIII|XIV|XV|XVI|MT)$'

printf "chrom\tstart\tend\tn_canonical_query_chromosomes\tn_other_canonical_query_chromosomes\n"

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT
windows_bed="$tmpdir/windows.bed"
hits_bedpe="$tmpdir/hits.bedpe"

for ((start = 0; start < length; start += window)); do
    end=$((start + window))
    if ((end > length)); then
        end=$length
    fi
    printf "%s\t%d\t%d\t%s:%d-%d\n" "$chrom" "$start" "$end" "$chrom" "$start" "$end" >> "$windows_bed"
done

impg query \
    -a "$syng" \
    --sequence-files "$agc" \
    -b "$windows_bed" \
    -d "$distance" \
    -o bedpe 2>/dev/null \
> "$hits_bedpe"

awk -v ref_chrom="$ref_chrom" -v chrom_re="$canonical_chrom_re" -v min_hit_len="$min_hit_len" '
    BEGIN { OFS = "\t" }
    FNR == NR {
        key = $4
        chrom[key] = $1
        start[key] = $2
        end[key] = $3
        order[++n] = key
        next
    }
    ($3 - $2) >= min_hit_len {
        query_name = $7
        split($1, q, "#")
        query_chrom = q[3]
        if (query_chrom ~ chrom_re) {
            seen[query_name SUBSEP query_chrom] = 1
            if (query_chrom != ref_chrom) {
                other[query_name SUBSEP query_chrom] = 1
            }
        }
    }
    END {
        for (k in seen) {
            split(k, p, SUBSEP)
            count[p[1]]++
        }
        for (k in other) {
            split(k, p, SUBSEP)
            other_count[p[1]]++
        }
        for (i = 1; i <= n; i++) {
            key = order[i]
            print chrom[key], start[key], end[key], count[key] + 0, other_count[key] + 0
        }
    }
' "$windows_bed" "$hits_bedpe"
