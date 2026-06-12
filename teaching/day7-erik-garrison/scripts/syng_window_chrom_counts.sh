#!/usr/bin/env bash
set -euo pipefail

chrom="${1:-SGDref#0#chrV}"
window="${2:-10000}"
distance="${3:-1k}"
syng="${4:-yeast235.syng}"
agc="${5:-yeast235.agc}"
fai="${6:-yeast235.fa.gz.fai}"

length=$(awk -v chrom="$chrom" '$1 == chrom {print $2}' "$fai")
if [[ -z "${length}" ]]; then
    echo "Could not find ${chrom} in ${fai}" >&2
    exit 1
fi

ref_chrom="${chrom##*#}"
canonical_chrom_re='^chr(I|II|III|IV|V|VI|VII|VIII|IX|X|XI|XII|XIII|XIV|XV|XVI|MT)$'

printf "chrom\tstart\tend\tn_canonical_query_chromosomes\tn_other_canonical_query_chromosomes\n"

for ((start = 0; start < length; start += window)); do
    end=$((start + window))
    if ((end > length)); then
        end=$length
    fi

    counts=$(
        impg query \
            -a "$syng" \
            --sequence-files "$agc" \
            -r "${chrom}:${start}-${end}" \
            -d "$distance" \
            -o bed 2>/dev/null \
        | cut -f 1 \
        | cut -f 3 -d '#' \
        | awk -v ref_chrom="$ref_chrom" -v chrom_re="$canonical_chrom_re" '
            $0 ~ chrom_re {
                seen[$0] = 1
                if ($0 != ref_chrom) {
                    other[$0] = 1
                }
            }
            END {
                for (c in seen) n++
                for (c in other) o++
                printf "%d\t%d", n + 0, o + 0
            }'
    )

    printf "%s\t%d\t%d\t%s\n" "$chrom" "$start" "$end" "$counts"
done
