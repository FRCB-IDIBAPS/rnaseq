#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# featureCounts QC: load *.featureCounts.tsv, build count matrix,
# and run a battery of QC checks (PCA + library metrics + diagnostics).
#
# Usage:
#   Rscript featurecounts_report.R [input_dir] [output_dir]
# Defaults: input_dir = ".", output_dir = "./qc_pca"
#
# Sample labels are clustered purely by string similarity (no design/
# metadata file) to keep same-looking replicates colour-consistent
# across plots -- this is a cosmetic grouping for readability, not a
# stand-in for an experimental design/condition assignment.
#
# Outputs (in output_dir):
#   merged_counts.tsv       gene x sample count matrix
#   qc_metrics.tsv          per-sample summary table
#   pca_scores.tsv          PC1..PC5 scores per sample
#   qc_report.pdf           ALL plots, one per page (combined report)
#   *.png                   one PNG per individual plot:
#     library_size, detected_genes, count_distributions,
#     sample_correlation, top_gene_fraction, assignment_rates,
#     pca_pc1_pc2, pca_pc3_pc4, scree_plot
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir  <- if (length(args) >= 1) args[1] else "."
output_dir <- if (length(args) >= 2) args[2] else file.path(input_dir, "qc_pca")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

TOP_N_GENES <- 50L   # for top-gene-fraction metric
CPM_DETECT  <- 1     # CPM threshold for "detected"

# Collect every plot here (in creation order) for the combined PDF report.
plot_list <- list()

# Write an individual PNG and register the plot for the combined PDF.
save_plot <- function(p, name, w = 7, h = 6, dpi = 150) {
  ggsave(file.path(output_dir, paste0(name, ".png")), p, width = w, height = h, dpi = dpi)
  plot_list[[name]] <<- p
  invisible(p)
}

# --- 1. Find files -------------------------------------------------
files <- list.files(input_dir, pattern = "\\.featureCounts\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No *.featureCounts.tsv files found in ", input_dir)
samples <- sub("\\.featureCounts\\.tsv$", "", basename(files))
message("Found ", length(files), " samples: ", paste(samples, collapse = ", "))

# --- 2. Load & merge into a gene x sample count matrix -------------
# nf-core/rnaseq featureCounts output: comment line starting with '#',
# header row, then columns: Geneid Chr Start End Strand Length <count>
read_fc <- function(f, s) {
  dt <- fread(f, skip = "Geneid")
  cnt_col <- ncol(dt)                       # count is the last column
  out <- dt[, .(Geneid, count = get(names(dt)[cnt_col]))]
  setnames(out, "count", s)
  out
}

mat_list <- Map(read_fc, files, samples)
counts   <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), mat_list)

gene_ids <- counts$Geneid
counts   <- as.matrix(counts[, -1])
rownames(counts) <- gene_ids
counts[is.na(counts)] <- 0
storage.mode(counts) <- "numeric"

fwrite(data.table(Geneid = rownames(counts), counts),
       file.path(output_dir, "merged_counts.tsv"), sep = "\t")
message("Merged matrix: ", nrow(counts), " genes x ", ncol(counts), " samples")

# --- Cosmetic sample grouping, purely from label string similarity -
# No design/metadata file is used here -- samples are clustered by
# normalised Levenshtein distance between their labels (e.g.
# "Ctrl_1"/"Ctrl_2"/"Ctrl-3" land in the same cluster with "Treated_1"
# in another) so that same-looking replicates share a plot colour.
# This is a readability aid only, never a biological condition call.
STEM_DIST_THRESHOLD <- 0.35  # normalised Levenshtein distance cut height

infer_label_clusters <- function(labels, threshold = STEM_DIST_THRESHOLD) {
  n <- length(labels)
  if (n <= 1) return(rep("Group1", n))

  d <- adist(labels)
  dimnames(d) <- list(labels, labels)
  norm <- outer(nchar(labels), nchar(labels), pmax)
  norm[norm == 0] <- 1L
  d_norm <- d / norm

  # "complete" linkage, plus forcing merge heights non-decreasing: "average"
  # linkage (and float error even with "complete") can otherwise produce a
  # non-monotonic dendrogram that makes cutree() error out on some label sets.
  hc <- hclust(as.dist(d_norm), method = "complete")
  hc$height <- cummax(hc$height)
  raw <- cutree(hc, h = threshold)          # aligned with `labels` order

  # Relabel clusters "Group1", "Group2", ... in order of first appearance
  # so labels stay stable/readable instead of hclust's internal numbering.
  ordered_ids <- unique(raw)
  paste0("Group", match(raw, ordered_ids))
}

group        <- infer_label_clusters(samples)
group_lookup <- setNames(group, samples)     # sample label -> cluster, for re-use below
group_levels <- unique(group)                # first-appearance order (not alphabetical --
                                              # "Group2"/"Group10" would sort wrong otherwise)

# Colour palette auto-scales to however many clusters were found, instead
# of a small hardcoded/dataset-specific colour vector.
qualitative_palette <- function(n) {
  pal_max <- brewer.pal.info["Set2", "maxcolors"]
  if (n <= 2) return(brewer.pal(3, "Set2")[seq_len(n)])   # brewer.pal needs n >= 3
  if (n <= pal_max) return(brewer.pal(n, "Set2"))
  colorRampPalette(brewer.pal(pal_max, "Set2"))(n)
}
group_colors <- setNames(qualitative_palette(length(group_levels)), group_levels)
grp_fill     <- scale_fill_manual(values = group_colors, name = "Label group")
grp_colour   <- scale_colour_manual(values = group_colors, name = "Label group")

# --- 3. Normalize (CPM + log2) ------------------------------------
keep     <- rowSums(counts) > 0
counts_f <- counts[keep, , drop = FALSE]
lib_size <- colSums(counts)                 # total assigned reads (all genes)
cpm      <- t(t(counts_f) / colSums(counts_f)) * 1e6
logcpm   <- log2(cpm + 1)

# Lightly filtered matrix for the distribution panel only: keep genes with
# CPM > CPM_DETECT in at least MIN_SAMPLES samples (= smallest group size).
# This lifts the density curves off the zero-inflation spike without
# affecting the other QC panels, which still see the full gene set.
MIN_SAMPLES  <- min(table(group))
expr_keep    <- rowSums(cpm > CPM_DETECT) >= MIN_SAMPLES
logcpm_expr  <- logcpm[expr_keep, , drop = FALSE]
message("Expression filter for distribution panel: ", sum(expr_keep),
        " / ", nrow(logcpm), " genes (CPM > ", CPM_DETECT,
        " in >= ", MIN_SAMPLES, " samples)")

# ==================================================================
# QC A. Library size (total assigned reads per sample)
# ==================================================================
lib_dt <- data.table(sample = samples, group = group, lib_size = lib_size)
p_lib <- ggplot(lib_dt, aes(reorder(sample, lib_size), lib_size / 1e6, fill = group)) +
  geom_col() +
  grp_fill +
  coord_flip() +
  labs(title = "Library size (assigned reads)", x = NULL,
       y = "Million assigned reads", fill = "Label group") +
  theme_bw()
save_plot(p_lib, "library_size")

# ==================================================================
# QC B. Detected genes per sample (CPM > threshold)
# ==================================================================
detected <- colSums(cpm > CPM_DETECT)
det_dt <- data.table(sample = samples, group = group, detected = detected)
p_det <- ggplot(det_dt, aes(reorder(sample, detected), detected, fill = group)) +
  geom_col() +
  grp_fill +
  coord_flip() +
  labs(title = sprintf("Detected genes (CPM > %g)", CPM_DETECT),
       x = NULL, y = "Number of genes", fill = "Label group") +
  theme_bw()
save_plot(p_det, "detected_genes")

# ==================================================================
# QC C. Per-sample count distributions (log2 CPM density)
# ==================================================================
long <- melt(data.table(logcpm_expr), measure.vars = colnames(logcpm_expr),
             variable.name = "sample", value.name = "logcpm")
long[, group := group_lookup[as.character(sample)]]
p_dist <- ggplot(long, aes(logcpm, colour = group, group = sample)) +
  geom_density(linewidth = 0.6, alpha = 0.8) +
  grp_colour +
  labs(title = "Count distributions (log2 CPM, expression-filtered)",
       x = "log2(CPM + 1)", y = "Density", colour = "Label group") +
  theme_bw()
save_plot(p_dist, "count_distributions")

# ==================================================================
# QC D. Sample-sample correlation heatmap (Spearman, clustered)
# ==================================================================
cor_mat <- cor(logcpm, method = "spearman")
hc      <- hclust(as.dist(1 - cor_mat))
lv      <- rownames(cor_mat)[hc$order]
cor_long <- melt(data.table(s1 = rownames(cor_mat), cor_mat),
                 id.vars = "s1", variable.name = "s2", value.name = "rho")
cor_long[, s1 := factor(s1, levels = lv)]
cor_long[, s2 := factor(s2, levels = lv)]
p_cor <- ggplot(cor_long, aes(s1, s2, fill = rho)) +
  geom_tile() +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b", name = "Spearman") +
  labs(title = "Sample-sample correlation (clustered)", x = NULL, y = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_cor, "sample_correlation")

# ==================================================================
# QC E. Top-gene fraction (% library in top N genes per sample)
# ==================================================================
top_frac <- apply(counts_f, 2, function(v) {
  sum(sort(v, decreasing = TRUE)[seq_len(min(TOP_N_GENES, length(v)))]) / sum(v) * 100
})
tf_dt <- data.table(sample = samples, group = group, top_frac = top_frac)
p_tf <- ggplot(tf_dt, aes(reorder(sample, top_frac), top_frac, fill = group)) +
  geom_col() +
  grp_fill +
  coord_flip() +
  labs(title = sprintf("Fraction of library in top %d genes", TOP_N_GENES),
       x = NULL, y = "% of assigned reads", fill = "Label group") +
  theme_bw()
save_plot(p_tf, "top_gene_fraction")

# ==================================================================
# QC F. featureCounts assignment rates (from .summary files)
# ==================================================================
summ_files <- paste0(files, ".summary")   # e.g. SAMPLE.featureCounts.tsv.summary
have_summ  <- file.exists(summ_files)
assign_rate <- rep(NA_real_, length(samples))
if (any(have_summ)) {
  summ_list <- list()
  for (i in which(have_summ)) {
    s  <- samples[i]
    dt <- fread(summ_files[i], header = TRUE)
    setnames(dt, c("Status", "Count"))
    dt[, sample := s][, group := group[i]]
    summ_list[[s]] <- dt
    tot <- dt[, sum(Count)]
    if (tot > 0) assign_rate[i] <- 100 * dt[Status == "Assigned", Count] / tot
  }
  summ_all <- rbindlist(summ_list)
  summ_all <- summ_all[, .(Count = sum(Count)), by = .(Status, sample, group)]
  summ_all <- summ_all[Status %in% summ_all[, sum(Count), by = Status][V1 > 0, Status]]
  p_asg <- ggplot(summ_all, aes(sample, Count / 1e6, fill = Status)) +
    geom_col() +
    coord_flip() +
    labs(title = "featureCounts read assignment", x = NULL,
         y = "Million reads", fill = "Status") +
    theme_bw()
  save_plot(p_asg, "assignment_rates")
} else {
  message("No .summary files found; skipping assignment-rate plot.")
}

# ==================================================================
# QC G. PCA (top variable genes) + scree + PC3/PC4
# ==================================================================
ntop <- min(10000L, nrow(logcpm))
vars <- apply(logcpm, 1, var)
top  <- order(vars, decreasing = TRUE)[seq_len(ntop)]
x    <- t(logcpm[top, , drop = FALSE])      # samples x genes

pca      <- prcomp(x, center = TRUE, scale. = FALSE)
var_expl <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
npc      <- min(5L, ncol(pca$x))

scores <- data.table(sample = rownames(pca$x), group = group)
for (k in seq_len(npc)) scores[[paste0("PC", k)]] <- pca$x[, k]
fwrite(scores, file.path(output_dir, "pca_scores.tsv"), sep = "\t")

pca_gg <- function(a, b) {
  ggplot(scores, aes(.data[[paste0("PC", a)]], .data[[paste0("PC", b)]],
                     colour = group, label = sample)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
    grp_colour +
    labs(x = sprintf("PC%d (%.1f%%)", a, var_expl[a]),
         y = sprintf("PC%d (%.1f%%)", b, var_expl[b]), colour = "Label group") +
    theme_bw()
}

p_pca12 <- pca_gg(1, 2) +
  ggtitle(sprintf("PCA QC (top %d variable genes, log2 CPM)", ntop))
save_plot(p_pca12, "pca_pc1_pc2")

# PC3/PC4 as its own PNG if enough components exist
if (npc >= 4) {
  save_plot(pca_gg(3, 4) + ggtitle("PCA QC (PC3 vs PC4)"), "pca_pc3_pc4")
}

# Scree plot
n_show  <- min(10L, length(var_expl))
scree   <- data.table(PC = factor(seq_len(n_show), levels = seq_len(n_show)),
                      var = var_expl[seq_len(n_show)])
p_scree <- ggplot(scree, aes(PC, var)) +
  geom_col(fill = "grey40") +
  geom_text(aes(label = sprintf("%.1f", var)), vjust = -0.4, size = 3) +
  labs(title = "Scree plot", x = "Principal component", y = "% variance explained") +
  theme_bw()
save_plot(p_scree, "scree_plot")

# ==================================================================
# Per-sample QC metrics table
# ==================================================================
qc_metrics <- data.table(
  sample          = samples,
  group           = group,
  lib_size        = lib_size,
  detected_genes  = detected,
  top50_pct       = round(top_frac, 2),
  assign_rate_pct = round(assign_rate, 2)
)
fwrite(qc_metrics, file.path(output_dir, "qc_metrics.tsv"), sep = "\t")

# ==================================================================
# Combined report: every plot, one per page, in a single PDF
# ==================================================================
report_pdf <- file.path(output_dir, "qc_report.pdf")
pdf(report_pdf, width = 8, height = 7, onefile = TRUE)
for (p in plot_list) print(p)
invisible(dev.off())
message("Combined report: ", report_pdf, " (", length(plot_list), " pages)")

message("Done. Outputs written to ", normalizePath(output_dir))
