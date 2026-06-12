BEGIN {
    OFS = "\t"
    if (window == "") {
        window = 10000
    }
    print "chrom", "start", "end", "n_canonical_query_chromosomes", "n_other_canonical_query_chromosomes"
}

function canonical_chrom(label) {
    return (label ~ /^chr(I|II|III|IV|V|VI|VII|VIII|IX|X|XI|XII|XIII|XIV|XV|XVI|MT)$/)
}

$6 ~ /^SGDref#0#/ {
    ref = $6
    ref_len = $7
    ref_start = $8
    ref_end = $9
    query = $1

    split(query, q, "#")
    query_chrom = q[3]
    if (query_chrom == "") {
        query_chrom = query
    }
    if (!canonical_chrom(query_chrom)) {
        next
    }
    split(ref, r, "#")
    ref_chrom = r[3]

    if (only_chrom != "" && ref != only_chrom) {
        next
    }

    first_bin = int(ref_start / window)
    last_bin = int((ref_end - 1) / window)
    for (bin = first_bin; bin <= last_bin; bin++) {
        key = ref SUBSEP bin
        seen[key SUBSEP query_chrom] = 1
        if (query_chrom != ref_chrom) {
            seen_other[key SUBSEP query_chrom] = 1
        }
        lengths[ref] = ref_len
    }
}

END {
    for (k in seen) {
        split(k, parts, SUBSEP)
        counts[parts[1] SUBSEP parts[2]]++
    }
    for (k in seen_other) {
        split(k, parts, SUBSEP)
        other_counts[parts[1] SUBSEP parts[2]]++
    }

    for (ref in lengths) {
        n_bins = int((lengths[ref] + window - 1) / window)
        for (bin = 0; bin < n_bins; bin++) {
            key = ref SUBSEP bin
            start = bin * window
            end = start + window
            if (end > lengths[ref]) {
                end = lengths[ref]
            }
            count = counts[key] + 0
            other_count = other_counts[key] + 0
            print ref, start, end, count, other_count
        }
    }
}
