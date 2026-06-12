args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: Rscript scripts/plot_window_counts.R counts.tsv output.png [title]", call. = FALSE)
}

input <- args[[1]]
output <- args[[2]]
plot_title <- if (length(args) >= 3) args[[3]] else "Query chromosomes per reference window"

suppressPackageStartupMessages({
  library(ggplot2)
})

x <- read.delim(input, check.names = FALSE)
x$midpoint_kb <- ((x$start + x$end) / 2) / 1000
x$chrom <- factor(x$chrom, levels = unique(x$chrom))
y_col <- if ("n_other_canonical_query_chromosomes" %in% names(x)) {
  "n_other_canonical_query_chromosomes"
} else {
  "n_query_chromosomes"
}
y_label <- if (y_col == "n_other_canonical_query_chromosomes") {
  "other canonical query chromosome labels"
} else {
  "distinct query chromosome labels"
}

p <- ggplot(x, aes(midpoint_kb, .data[[y_col]])) +
  geom_col(width = 9, fill = "#3b6fb6") +
  facet_wrap(~ chrom, scales = "free_x", ncol = 1) +
  labs(
    title = plot_title,
    x = "SGDref position (kb)",
    y = y_label
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    panel.grid.minor = element_blank()
  )

ggsave(output, p, width = 10, height = max(3, 1.5 * length(unique(x$chrom))), dpi = 160, limitsize = FALSE)
