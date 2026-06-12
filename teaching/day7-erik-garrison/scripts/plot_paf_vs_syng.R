args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: Rscript scripts/plot_paf_vs_syng.R counts.tsv output.png [title]", call. = FALSE)
}

input <- args[[1]]
output <- args[[2]]
plot_title <- if (length(args) >= 3) args[[3]] else "PAF vs syng window counts"

suppressPackageStartupMessages({
  library(ggplot2)
})

x <- read.delim(input, check.names = FALSE)
x$midpoint_kb <- ((x$start + x$end) / 2) / 1000
y_field <- if ("paf_other_count" %in% names(x) && "syng_other_count" %in% names(x)) {
  c("paf_other_count", "syng_other_count")
} else {
  c("paf_count", "syng_count")
}
y_label <- if (identical(y_field, c("paf_other_count", "syng_other_count"))) {
  "other canonical query chromosome labels"
} else {
  "distinct query chromosome labels"
}

long <- rbind(
  data.frame(chrom = x$chrom, midpoint_kb = x$midpoint_kb, backend = "PAF", n_query_chromosomes = x[[y_field[[1]]]]),
  data.frame(chrom = x$chrom, midpoint_kb = x$midpoint_kb, backend = "syng", n_query_chromosomes = x[[y_field[[2]]]])
)

p <- ggplot(long, aes(midpoint_kb, n_query_chromosomes, color = backend)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  facet_wrap(~ chrom, scales = "free_x", ncol = 1) +
  scale_color_manual(values = c(PAF = "#3b6fb6", syng = "#b64a3b")) +
  labs(
    title = plot_title,
    x = "SGDref position (kb)",
    y = y_label,
    color = "backend"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(output, p, width = 10, height = 3.5, dpi = 160, limitsize = FALSE)
