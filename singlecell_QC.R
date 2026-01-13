###############################################################################
#QC-only 스크립트 
# [1] 300 ≤ nFeature ≤ 6000
# [2] Mito ≤ 10%
# [3] 하위 2% 제거
# [4] Complexity(log10GenesPerUMI) 하위 2% 제거: 0.80–0.90 clamp0.80–0.90 clamp
###############################################################################


#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

options(Seurat.object.assay.version = "v5")

# =========================
# Config (QC-only)
# =========================
IN_RDS  <- "/data1/pubdata/ipf_scrna-seq/ipf_obj_merged.rds"
OUT_RDS <- "/data1/pubdata/ipf_scrna-seq/ipf_obj_merged.qc_filtered.QConly.rds"
OUT_QC_SUMMARY <- "/data1/pubdata/ipf_scrna-seq/ipf_obj_merged.qc_summary.QConly.tsv"
PLOT_DIR <- "/data1/pubdata/ipf_scrna-seq/qc_plots_compare_QConly"

# QC cutoffs
NF_LO <- 300
NF_HI <- 6000
MT_HI <- 10  # percent

# Complexity: log10GenesPerUMI = log10(nFeature+1)/log10(nCount+1)
LOWQ_P   <- 0.02     # bottom 2% quantile
CLAMP_LO <- 0.80
CLAMP_HI <- 0.90

# =========================
# Utils
# =========================
clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

safe_name <- function(x) {
  x <- gsub("[/\\s]+", "_", x)
  gsub("[^A-Za-z0-9._-]", "_", x)
}
dir_create <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# Prevent "cannot xtfrm data frames" by flattening meta columns to atomic vectors
fix_meta_atomic <- function(seu) {
  md <- as.data.frame(seu@meta.data, stringsAsFactors = FALSE)
  for (nm in colnames(md)) {
    x <- md[[nm]]
    if (is.data.frame(x)) {
      md[[nm]] <- x[[1]]
      next
    }
    if (is.matrix(x) || is.array(x)) {
      md[[nm]] <- as.vector(x)
      next
    }
    if (is.list(x)) {
      if (all(vapply(x, function(z) length(z) == 1, logical(1)))) {
        md[[nm]] <- unlist(x)
      } else {
        md[[nm]] <- vapply(x, function(z) paste(z, collapse=";"), character(1))
      }
      next
    }
  }
  colnames(md) <- make.unique(colnames(md))
  seu@meta.data <- md
  seu
}

# Remove any DoubletFinder leftovers (even though QC-only, safer for re-runs)
drop_df_cols <- function(seu) {
  md <- seu@meta.data
  drop <- grep("^(pANN_|DF\\.classifications|DF\\.pANN|DF\\.|DF_class$)", colnames(md), value = TRUE)
  if (length(drop) > 0) md <- md[, setdiff(colnames(md), drop), drop = FALSE]
  if (any(duplicated(colnames(md)))) colnames(md) <- make.unique(colnames(md))
  seu@meta.data <- md
  seu
}

pct <- function(n, denom) ifelse(denom == 0, NA_real_, n / denom * 100)

primary_driver_label <- function(only_nf, only_mt, only_cx) {
  v <- c(nFeature = only_nf, mt = only_mt, complexity = only_cx)
  nm <- names(v)[which.max(v)]
  if (sum(v == max(v)) > 1) return("tie")
  nm
}

# =========================
# QC metrics
# =========================
ensure_percent_mt <- function(obj) {
  if (!"percent.mt" %in% colnames(obj@meta.data)) {
    feats <- rownames(obj)
    if (any(grepl("^MT-", feats))) {
      obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    } else if (any(grepl("^mt-", feats))) {
      obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
    } else {
      warning("No MT- or mt- genes detected; setting percent.mt = 0")
      obj$percent.mt <- 0
    }
  }
  obj
}

add_complexity <- function(obj) {
  if (!"log10GenesPerUMI" %in% colnames(obj@meta.data)) {
    obj$log10GenesPerUMI <- log10(obj$nFeature_RNA + 1) / log10(obj$nCount_RNA + 1)
  }
  obj
}

compute_complexity_cutoff <- function(obj) {
  x <- obj$log10GenesPerUMI
  q <- as.numeric(stats::quantile(x, probs = LOWQ_P, na.rm = TRUE))
  clamp(q, CLAMP_LO, CLAMP_HI)
}

# =========================
# Plot helpers (violin + hist+density + scatter, before/after)
# =========================
qc_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey30"),
      axis.title = element_text(face = "bold")
    )
}
COL_FILL   <- "#4C78A8"
COL_DENS   <- "#F58518"
COL_CUTOFF <- "#E45756"
COL_POINT  <- "#54A24B"
COL_FAIL   <- "#E45756"

violin_plot <- function(df, metric, cut_lines = numeric(), title = NULL, subtitle = NULL) {
  p <- ggplot(df, aes(x = "all", y = .data[[metric]])) +
    geom_violin(fill = COL_FILL, color = NA, alpha = 0.85, trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.size = 0.2, fill = "white", color = "grey20") +
    labs(x = NULL, y = metric, title = title, subtitle = subtitle) +
    qc_theme() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  for (cl in cut_lines) {
    p <- p + geom_hline(yintercept = cl, linetype = "dashed", color = COL_CUTOFF, linewidth = 0.8)
  }
  p
}

hist_plot <- function(df, metric, cut_lines = numeric(), title = NULL, subtitle = NULL, bins = 60) {
  p <- ggplot(df, aes(x = .data[[metric]])) +
    geom_histogram(aes(y = after_stat(density)), bins = bins, fill = COL_FILL, alpha = 0.65, color = "white") +
    geom_density(color = COL_DENS, linewidth = 1.0, adjust = 1.1) +
    labs(x = metric, y = "Density", title = title, subtitle = subtitle) +
    qc_theme()
  for (cl in cut_lines) {
    p <- p + geom_vline(xintercept = cl, linetype = "dashed", color = COL_CUTOFF, linewidth = 0.8)
  }
  p
}

scatter_plot_passfail <- function(df, x, y, vlines = numeric(), hlines = numeric(),
                                  title = NULL, subtitle = NULL,
                                  alpha = 0.45, size = 0.55) {
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], color = pass_qc)) +
    geom_point(alpha = alpha, size = size) +
    scale_color_manual(values = c(`TRUE` = COL_POINT, `FALSE` = COL_FAIL)) +
    labs(x = x, y = y, title = title, subtitle = subtitle, color = "Pass QC") +
    qc_theme()
  for (cl in vlines) p <- p + geom_vline(xintercept = cl, linetype="dashed", color=COL_CUTOFF, linewidth=0.8)
  for (cl in hlines) p <- p + geom_hline(yintercept = cl, linetype="dashed", color=COL_CUTOFF, linewidth=0.8)
  p
}

save_side_by_side <- function(p_left, p_right, out_png, width = 13.8, height = 5.2, dpi = 180) {
  ggsave(out_png, p_left + p_right, width = width, height = height, dpi = dpi)
}

make_compare_plots <- function(obj_before, obj_after, sample_id, out_dir, cx_cut,
                               NF_LO, NF_HI, MT_HI) {
  df_before <- obj_before@meta.data %>%
    dplyr::transmute(
      nFeature_RNA = nFeature_RNA,
      nCount_RNA = nCount_RNA,
      percent.mt = percent.mt,
      log10GenesPerUMI = log10GenesPerUMI,
      pass_qc = (nFeature_RNA >= NF_LO) & (nFeature_RNA <= NF_HI) &
        (percent.mt <= MT_HI) &
        (log10GenesPerUMI >= cx_cut)
    )
  df_after <- obj_after@meta.data %>%
    dplyr::transmute(
      nFeature_RNA = nFeature_RNA,
      nCount_RNA = nCount_RNA,
      percent.mt = percent.mt,
      log10GenesPerUMI = log10GenesPerUMI,
      pass_qc = TRUE
    )
  
  sdir <- file.path(out_dir, safe_name(sample_id))
  dir_create(sdir)
  
  # Violin
  save_side_by_side(
    violin_plot(df_before, "nFeature_RNA", c(NF_LO, NF_HI),
                title = paste0(sample_id, " | nFeature_RNA"), subtitle = "Before"),
    violin_plot(df_after,  "nFeature_RNA", c(NF_LO, NF_HI),
                title = paste0(sample_id, " | nFeature_RNA"), subtitle = "After"),
    file.path(sdir, "compare_violin_nFeature_RNA.png"),
    height = 4.8
  )
  
  save_side_by_side(
    violin_plot(df_before, "percent.mt", c(MT_HI),
                title = paste0(sample_id, " | percent.mt"), subtitle = "Before"),
    violin_plot(df_after,  "percent.mt", c(MT_HI),
                title = paste0(sample_id, " | percent.mt"), subtitle = "After"),
    file.path(sdir, "compare_violin_percent.mt.png"),
    height = 4.8
  )
  
  save_side_by_side(
    violin_plot(df_before, "log10GenesPerUMI", c(cx_cut),
                title = paste0(sample_id, " | log10GenesPerUMI"),
                subtitle = sprintf("Before (cutoff=%.3f)", cx_cut)),
    violin_plot(df_after,  "log10GenesPerUMI", c(cx_cut),
                title = paste0(sample_id, " | log10GenesPerUMI"), subtitle = "After"),
    file.path(sdir, "compare_violin_log10GenesPerUMI.png"),
    height = 4.8
  )
  
  # Hist + density
  save_side_by_side(
    hist_plot(df_before, "nFeature_RNA", c(NF_LO, NF_HI),
              title = paste0(sample_id, " | nFeature_RNA (hist+density)"), subtitle = "Before"),
    hist_plot(df_after,  "nFeature_RNA", c(NF_LO, NF_HI),
              title = paste0(sample_id, " | nFeature_RNA (hist+density)"), subtitle = "After"),
    file.path(sdir, "compare_hist_nFeature_RNA.png"),
    height = 4.8
  )
  
  save_side_by_side(
    hist_plot(df_before, "percent.mt", c(MT_HI),
              title = paste0(sample_id, " | percent.mt (hist+density)"), subtitle = "Before"),
    hist_plot(df_after,  "percent.mt", c(MT_HI),
              title = paste0(sample_id, " | percent.mt (hist+density)"), subtitle = "After"),
    file.path(sdir, "compare_hist_percent.mt.png"),
    height = 4.8
  )
  
  save_side_by_side(
    hist_plot(df_before, "log10GenesPerUMI", c(cx_cut),
              title = paste0(sample_id, " | log10GenesPerUMI (hist+density)"),
              subtitle = sprintf("Before (cutoff=%.3f)", cx_cut)),
    hist_plot(df_after,  "log10GenesPerUMI", c(cx_cut),
              title = paste0(sample_id, " | log10GenesPerUMI (hist+density)"), subtitle = "After"),
    file.path(sdir, "compare_hist_log10GenesPerUMI.png"),
    height = 4.8
  )
  
  # Scatter
  save_side_by_side(
    scatter_plot_passfail(df_before, "nCount_RNA", "nFeature_RNA",
                          hlines = c(NF_LO, NF_HI),
                          title = paste0(sample_id, " | nCount vs nFeature"),
                          subtitle = "Before"),
    scatter_plot_passfail(df_after, "nCount_RNA", "nFeature_RNA",
                          hlines = c(NF_LO, NF_HI),
                          title = paste0(sample_id, " | nCount vs nFeature"),
                          subtitle = "After"),
    file.path(sdir, "compare_scatter_nCount_vs_nFeature.png"),
    height = 5.2
  )
}

# =========================
# MAIN
# =========================
dir_create(PLOT_DIR)
cat(sprintf("[MAIN] reading: %s\n", IN_RDS))
obj <- readRDS(IN_RDS)

# sample column auto-detect
candidate_cols <- c("Sample", "orig.ident", "sample", "SampleName", "SampleID")
sample_col <- candidate_cols[candidate_cols %in% colnames(obj@meta.data)][1]
if (is.na(sample_col) || length(sample_col) == 0) {
  stop(sprintf(
    "No sample column found in obj@meta.data. Tried: %s\nAvailable columns (first 30): %s",
    paste(candidate_cols, collapse = ", "),
    paste(head(colnames(obj@meta.data), 30), collapse = ", ")
  ))
}
cat(sprintf("[MAIN] using sample column: %s\n", sample_col))

# sanitize meta
obj <- fix_meta_atomic(obj)
obj <- drop_df_cols(obj)

all_samples <- sort(unique(obj@meta.data[[sample_col]]))
cat(sprintf("[MAIN] total samples = %d\n", length(all_samples)))

summary_rows  <- list()
filtered_list <- list()

for (nm in all_samples) {
  cat("\n=============================================\n")
  cat(sprintf("=== Sample: %s ===\n", as.character(nm)))
  
  cells_use <- rownames(obj@meta.data)[obj@meta.data[[sample_col]] == nm]
  if (length(cells_use) == 0) {
    message(sprintf("[WARN] no cells for sample=%s -> skip", as.character(nm)))
    next
  }
  o <- subset(obj, cells = cells_use)
  
  o <- fix_meta_atomic(o)
  o <- drop_df_cols(o)
  o <- ensure_percent_mt(o)
  o <- add_complexity(o)
  
  cx_cut <- compute_complexity_cutoff(o)
  cat(sprintf("Complexity cutoff (%.0f%% quantile, clamped): %.4f\n", LOWQ_P*100, cx_cut))
  
  n_before <- ncol(o)
  
  # Per-cell fail flags (QC-only)
  fail_nf <- (o$nFeature_RNA < NF_LO) | (o$nFeature_RNA > NF_HI)
  fail_mt <- (o$percent.mt > MT_HI)
  fail_cx <- (o$log10GenesPerUMI < cx_cut)
  
  keep <- !(fail_nf | fail_mt | fail_cx)
  o_f <- subset(o, cells = colnames(o)[keep])
  
  n_after <- ncol(o_f)
  removed <- n_before - n_after
  removed_pct <- ifelse(n_before == 0, NA_real_, removed / n_before * 100)
  cat(sprintf("Cells: %d -> %d (removed %d, %.2f%%)\n", n_before, n_after, removed, removed_pct))
  
  # Breakdown counts
  fail_nFeature_n <- sum(fail_nf, na.rm = TRUE)
  fail_mt_n       <- sum(fail_mt, na.rm = TRUE)
  fail_complexity_n <- sum(fail_cx, na.rm = TRUE)
  
  only_nf <- sum(fail_nf & !fail_mt & !fail_cx, na.rm = TRUE)
  only_mt <- sum(!fail_nf & fail_mt & !fail_cx, na.rm = TRUE)
  only_cx <- sum(!fail_nf & !fail_mt & fail_cx, na.rm = TRUE)
  
  nf_mt <- sum(fail_nf & fail_mt & !fail_cx, na.rm = TRUE)
  nf_cx <- sum(fail_nf & !fail_mt & fail_cx, na.rm = TRUE)
  mt_cx <- sum(!fail_nf & fail_mt & fail_cx, na.rm = TRUE)
  all3  <- sum(fail_nf & fail_mt & fail_cx, na.rm = TRUE)
  
  # Driver
  primary_driver <- primary_driver_label(only_nf, only_mt, only_cx)
  exclusive_sum <- only_nf + only_mt + only_cx
  
  primary_driver_removed_pct <- {
    if (primary_driver == "tie") NA_real_ else
      pct(c(nFeature = only_nf, mt = only_mt, complexity = only_cx)[[primary_driver]], removed)
  }
  primary_driver_exclusive_pct <- {
    if (primary_driver == "tie") NA_real_ else
      pct(c(nFeature = only_nf, mt = only_mt, complexity = only_cx)[[primary_driver]], exclusive_sum)
  }
  
  # Shares (removed 기반)
  mt_only_share_removed_pct <- pct(only_mt, removed)
  nf_only_share_removed_pct <- pct(only_nf, removed)
  cx_only_share_removed_pct <- pct(only_cx, removed)
  
  nf_mt_share_removed_pct <- pct(nf_mt, removed)
  nf_cx_share_removed_pct <- pct(nf_cx, removed)
  mt_cx_share_removed_pct <- pct(mt_cx, removed)
  all3_share_removed_pct  <- pct(all3, removed)
  
  # plots
  tryCatch({
    make_compare_plots(
      obj_before = o,
      obj_after  = o_f,
      sample_id  = as.character(nm),
      out_dir    = PLOT_DIR,
      cx_cut     = cx_cut,
      NF_LO      = NF_LO,
      NF_HI      = NF_HI,
      MT_HI      = MT_HI
    )
  }, error = function(e) {
    message(sprintf("[PLOT] failed (continue): %s", conditionMessage(e)))
  })
  
  # summary
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Sample = as.character(nm),
    nCells_before = n_before,
    nCells_after  = n_after,
    removed       = removed,
    removed_pct   = removed_pct,
    
    nFeature_lo   = NF_LO,
    nFeature_hi   = NF_HI,
    percent_mt_hi = MT_HI,
    complexity_quantile = LOWQ_P,
    complexity_cutoff = cx_cut,
    
    fail_nFeature_n = fail_nFeature_n,
    fail_nFeature_pct = pct(fail_nFeature_n, n_before),
    fail_mt_n = fail_mt_n,
    fail_mt_pct = pct(fail_mt_n, n_before),
    fail_complexity_n = fail_complexity_n,
    fail_complexity_pct = pct(fail_complexity_n, n_before),
    
    fail_only_nFeature_n = only_nf,
    fail_only_mt_n = only_mt,
    fail_only_complexity_n = only_cx,
    fail_nFeature_and_mt_n = nf_mt,
    fail_nFeature_and_complexity_n = nf_cx,
    fail_mt_and_complexity_n = mt_cx,
    fail_all3_n = all3,
    
    primary_driver = primary_driver,
    primary_driver_removed_pct = primary_driver_removed_pct,
    primary_driver_exclusive_pct = primary_driver_exclusive_pct,
    
    mt_only_share_removed_pct = mt_only_share_removed_pct,
    nf_only_share_removed_pct = nf_only_share_removed_pct,
    cx_only_share_removed_pct = cx_only_share_removed_pct,
    
    nf_mt_share_removed_pct = nf_mt_share_removed_pct,
    nf_cx_share_removed_pct = nf_cx_share_removed_pct,
    mt_cx_share_removed_pct = mt_cx_share_removed_pct,
    all3_share_removed_pct  = all3_share_removed_pct,
    
    stringsAsFactors = FALSE
  )
  
  filtered_list[[as.character(nm)]] <- o_f
  
  rm(o, o_f, cells_use)
  gc()
}

# merge
cat("\n[MERGE] merging filtered objects...\n")
if (length(filtered_list) == 0) stop("No filtered objects to merge.")

nm0 <- names(filtered_list)[1]
obj_filtered <- filtered_list[[nm0]]
if (length(filtered_list) > 1) {
  for (nm in names(filtered_list)[-1]) {
    obj_filtered <- merge(obj_filtered, y = filtered_list[[nm]], merge.data = FALSE)
  }
}

# save
summary_df <- do.call(rbind, summary_rows)
write.table(summary_df, file = OUT_QC_SUMMARY, sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(obj_filtered, OUT_RDS)

cat(sprintf("[DONE] saved filtered object: %s\n", OUT_RDS))
cat(sprintf("[DONE] saved summary: %s\n", OUT_QC_SUMMARY))
cat(sprintf("[DONE] plots dir: %s\n", PLOT_DIR))





###############################################################################
# 1. removed_pct 기준 Top10 샘플 자동 추출 (TSV 저장)
# 2. Top10 샘플의 제거 원인 구성(단독/교집합) stacked bar plot (PNG 저장)
# 3. Top10 샘플의 주범(primary_driver) + 주범 비율(primary_driver_removed_pct) 표 (TSV 저장)
###############################################################################

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

# =========================
# User inputs
# =========================
QC_SUMMARY_TSV <- "/data1/pubdata/ipf_scrna-seq/ipf_obj_merged.qc_summary.QConly.tsv"   # e.g. "/data1/pubdata/ipf_scrna-seq/qc/qc_summary_final.tsv"
OUTDIR <- "/data1/pubdata/ipf_scrna-seq/"          # e.g. "/data1/pubdata/ipf_scrna-seq/qc/qc_summary_plots"
TOPN <- 10

# =========================
# Helpers
# =========================
stop_if_empty <- function(x, name) {
  if (is.null(x) || !nzchar(x)) stop(sprintf("Please set %s path.", name))
}

required_cols <- c(
  "Sample","nCells_before","nCells_after","removed","removed_pct",
  "primary_driver","primary_driver_removed_pct","primary_driver_exclusive_pct",
  "mt_only_share_removed_pct","nf_only_share_removed_pct","cx_only_share_removed_pct",
  "nf_mt_share_removed_pct","nf_cx_share_removed_pct","mt_cx_share_removed_pct","all3_share_removed_pct",
  "fail_only_mt_n","fail_only_nFeature_n","fail_only_complexity_n",
  "fail_nFeature_and_mt_n","fail_nFeature_and_complexity_n","fail_mt_and_complexity_n","fail_all3_n"
)

# =========================
# Main
# =========================
stop_if_empty(QC_SUMMARY_TSV, "QC_SUMMARY_TSV")
stop_if_empty(OUTDIR, "OUTDIR")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

df <- readr::read_tsv(QC_SUMMARY_TSV, show_col_types = FALSE)

missing <- setdiff(required_cols, colnames(df))
if (length(missing) > 0) {
  stop("QC summary TSV is missing required columns:\n  - ",
       paste(missing, collapse = "\n  - "))
}

# ensure numeric
num_cols <- c(
  "removed_pct",
  "primary_driver_removed_pct","primary_driver_exclusive_pct",
  "mt_only_share_removed_pct","nf_only_share_removed_pct","cx_only_share_removed_pct",
  "nf_mt_share_removed_pct","nf_cx_share_removed_pct","mt_cx_share_removed_pct","all3_share_removed_pct"
)
df <- df %>%
  mutate(across(all_of(num_cols), as.numeric))

# TopN by removed_pct
top <- df %>%
  arrange(desc(removed_pct), desc(removed)) %>%
  slice_head(n = TOPN)

# Save TopN table
top_path <- file.path(OUTDIR, sprintf("top%d_by_removed_pct.tsv", TOPN))
readr::write_tsv(top, top_path)

# Save driver-focused compact table
top_driver <- top %>%
  transmute(
    Sample,
    nCells_before, nCells_after,
    removed, removed_pct,
    primary_driver,
    primary_driver_removed_pct,
    primary_driver_exclusive_pct,
    mt_only_share_removed_pct,
    nf_only_share_removed_pct,
    cx_only_share_removed_pct,
    nf_mt_share_removed_pct,
    nf_cx_share_removed_pct,
    mt_cx_share_removed_pct,
    all3_share_removed_pct
  )
driver_path <- file.path(OUTDIR, sprintf("top%d_driver_breakdown.tsv", TOPN))
readr::write_tsv(top_driver, driver_path)

# ---------- Plot 1: stacked bar by removed composition (percent of removed)
# We'll stack these shares (already percent-of-removed) per sample:
#   mt_only, nf_only, cx_only, nf_mt, nf_cx, mt_cx, all3
comp_long <- top %>%
  select(
    Sample,
    mt_only_share_removed_pct,
    nf_only_share_removed_pct,
    cx_only_share_removed_pct,
    nf_mt_share_removed_pct,
    nf_cx_share_removed_pct,
    mt_cx_share_removed_pct,
    all3_share_removed_pct
  ) %>%
  pivot_longer(
    cols = -Sample,
    names_to = "component",
    values_to = "pct"
  ) %>%
  mutate(
    component = recode(component,
                       mt_only_share_removed_pct = "mt only",
                       nf_only_share_removed_pct = "nFeature only",
                       cx_only_share_removed_pct = "complexity only",
                       nf_mt_share_removed_pct = "nFeature + mt",
                       nf_cx_share_removed_pct = "nFeature + complexity",
                       mt_cx_share_removed_pct = "mt + complexity",
                       all3_share_removed_pct = "all 3"
    ),
    Sample = factor(Sample, levels = top$Sample)
  )

p1 <- ggplot(comp_long, aes(x = Sample, y = pct, fill = component)) +
  geom_col(width = 0.85) +
  coord_flip() +
  labs(
    title = sprintf("Top %d samples by removed_pct: composition of removed cells", TOPN),
    x = NULL,
    y = "% of removed cells",
    fill = "Fail overlap"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

p1_path <- file.path(OUTDIR, sprintf("top%d_removed_composition_stacked.png", TOPN))
ggsave(p1_path, p1, width = 12, height = 7, dpi = 200)

# ---------- Plot 2: removed_pct bar + label primary_driver (quick view)
p2 <- top %>%
  mutate(Sample = factor(Sample, levels = top$Sample)) %>%
  ggplot(aes(x = Sample, y = removed_pct, fill = primary_driver)) +
  geom_col(width = 0.85) +
  coord_flip() +
  geom_text(
    aes(label = sprintf("%s (%.1f%%)", primary_driver, primary_driver_removed_pct)),
    hjust = -0.05,
    size = 3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = sprintf("Top %d samples by removed_pct: removed_pct and primary driver", TOPN),
    x = NULL,
    y = "removed_pct",
    fill = "primary_driver"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

p2_path <- file.path(OUTDIR, sprintf("top%d_removed_pct_primary_driver.png", TOPN))
ggsave(p2_path, p2, width = 12, height = 6, dpi = 200)

# ---------- (Optional) Plot 3: counts-based overlap composition (raw n), not percent
# This is useful if you want absolute cell numbers in overlaps.
counts_long <- top %>%
  select(
    Sample,
    fail_only_mt_n,
    fail_only_nFeature_n,
    fail_only_complexity_n,
    fail_nFeature_and_mt_n,
    fail_nFeature_and_complexity_n,
    fail_mt_and_complexity_n,
    fail_all3_n
  ) %>%
  pivot_longer(
    cols = -Sample,
    names_to = "component",
    values_to = "n"
  ) %>%
  mutate(
    component = recode(component,
                       fail_only_mt_n = "mt only",
                       fail_only_nFeature_n = "nFeature only",
                       fail_only_complexity_n = "complexity only",
                       fail_nFeature_and_mt_n = "nFeature + mt",
                       fail_nFeature_and_complexity_n = "nFeature + complexity",
                       fail_mt_and_complexity_n = "mt + complexity",
                       fail_all3_n = "all 3"
    ),
    Sample = factor(Sample, levels = top$Sample)
  )

p3 <- ggplot(counts_long, aes(x = Sample, y = n, fill = component)) +
  geom_col(width = 0.85) +
  coord_flip() +
  labs(
    title = sprintf("Top %d samples by removed_pct: composition of removed cells (counts)", TOPN),
    x = NULL,
    y = "Number of removed cells",
    fill = "Fail overlap"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

p3_path <- file.path(OUTDIR, sprintf("top%d_removed_composition_counts_stacked.png", TOPN))
ggsave(p3_path, p3, width = 12, height = 7, dpi = 200)

message("[DONE]")
message("Top table: ", top_path)
message("Driver table: ", driver_path)
message("Plot1 (removed composition %): ", p1_path)
message("Plot2 (removed_pct + driver): ", p2_path)
message("Plot3 (removed composition counts): ", p3_path)






###############################################################################

# Visualize QC filtered sampee dist plot

###############################################################################
# 필요한 라이브러리 로드
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra) # 여러 그래프를 한 페이지에 배치하기 위함

# 1. 데이터 로드
file_path <- "/data1/pubdata/ipf_scrna-seq/ipf_obj_merged.qc_summary.QConly.tsv"
df <- read.table(file_path, header = TRUE, sep = "\t")

# 2. 시각화를 위한 데이터 정렬 (제거율 기준 내림차순)
df <- df %>% arrange(desc(removed_pct))
# Factor 레벨을 고정하여 그래프에서 정렬된 순서 유지
df$Sample <- factor(df$Sample, levels = df$Sample)

# --- 그래프 1: Stacked Bar Chart (세포 수 비교) ---
# 가독성을 위해 데이터를 'long' 포맷으로 변경
df_long <- df %>%
  select(Sample, nCells_after, removed) %>%
  pivot_longer(cols = c(nCells_after, removed), names_to = "Status", values_to = "Count")

# Status 레이블 변경 (범례용)
df_long$Status <- factor(df_long$Status,
                         levels = c("removed", "nCells_after"),
                         labels = c("Removed", "Passed"))

p1 <- ggplot(df_long, aes(x = Sample, y = Count, fill = Status)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Removed" = "#C44E52", "Passed" = "#4C72B0")) +
  theme_minimal() +
  labs(title = "Cell Count Comparison (Before vs After QC)",
       x = NULL, y = "Number of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- 그래프 2: Bar Chart (제거율 %) ---
p2 <- ggplot(df, aes(x = Sample, y = removed_pct)) +
  geom_bar(stat = "identity", fill = "#E67E22") +
  geom_hline(yintercept = mean(df$removed_pct), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "Percentage of Cells Removed per Sample",
       subtitle = paste0("Mean Removal: ", round(mean(df$removed_pct), 1), "%"),
       x = "Sample", y = "Removed Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3. 두 그래프를 수직으로 결합하여 출력 및 저장
combined_plot <- grid.arrange(p1, p2, ncol = 1)

# 파일로 저장
ggsave("qc_summary_plots_R.png", combined_plot, width = 12, height = 12, dpi = 300)


#############################################

# Drop QC Fail Sample

#############################################
library(Seurat)
library(dplyr)
library(readr)

# 1. 샘플 제외 체크표 로드
# 이 표의 Final_Status 컬럼을 기준으로 제외 대상을 결정합니다.
exclusion_table <- read_delim("sample_exclusion_check_table.tsv", delim = "\t")

# 2. 제외(Dropped)할 샘플 리스트 추출
dropped_samples <- exclusion_table %>%
  filter(Final_Status == "Dropped") %>%
  pull(Sample)

cat("제거될 샘플 리스트:", paste(dropped_samples, collapse = ", "), "\n")

# 3. 원본 RDS 파일 로드
# 파일명: ipf_obj_merged.qc_filtered.QConly.rds
rds_path <- "ipf_obj_merged.qc_filtered.QConly.rds"
seurat_obj <- readRDS(rds_path)

# 4. Dropped 샘플 제거 (Subset 활용)
# metadata의 'orig.ident' 또는 'Sample' 컬럼에 샘플명이 저장되어 있다고 가정합니다.
# (사용자의 metadata 컬럼명에 맞춰 수정이 필요할 수 있습니다. 여기서는 'Sample' 기준)
if ("Sample" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj_filtered <- subset(seurat_obj, subset = Sample %in% dropped_samples, invert = TRUE)
} else {
  # 만약 'Sample' 컬럼이 없다면 'orig.ident'를 기본으로 시도합니다.
  seurat_obj_filtered <- subset(seurat_obj, subset = orig.ident %in% dropped_samples, invert = TRUE)
}

# 5. 결과 확인
cat("필터링 전 세포 수:", ncol(seurat_obj), "\n")
cat("필터링 후 세포 수:", ncol(seurat_obj_filtered), "\n")
cat("남은 샘플 수:", length(unique(seurat_obj_filtered@meta.data$Sample)), "\n")

# 6. 새로운 RDS 파일로 저장
# 파일명: ipf_obj_merged.qc_filtered.QC_dropped.rds
saveRDS(seurat_obj_filtered, "ipf_obj_merged.qc_filtered.QC_Fail_dropped.rds")

cat("성공적으로 필터링된 객체가 'ipf_obj_merged.qc_filtered.QC_dropped.rds'로 저장되었습니다.\n")



#############################################

# Doubletfinder

#############################################
#!/usr/bin/env Rscript

# 메모리 제한 확장 (1024GiB 서버 사양 반영)
options(future.globals.maxSize = 50 * 1024^3) # 프로세스당 최대 50GB 허용

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(dplyr)
  library(future)
  library(future.apply)
})

# =========================
# 병렬 처리 설정
# =========================
# 서버 여유 자원을 고려하여 24개 프로세스로 설정
plan("multisession", workers = 24) 
cat(sprintf("[MAIN] Parallel mode enabled: %d workers\n", 24))

# =========================
# Config
# =========================
IN_RDS              <- "/data1/pubdata/ipf_scrna-seq/ipf_obj_merged.qc_filtered.QC_Fail_dropped.rds"
OUT_RDS_WITH_DF     <- "/data1/pubdata/ipf_scrna-seq/ipf_obj_merged.qc_filtered.QConly_with_DF.rds"
OUT_RDS_DF_FILTERED <- "/data1/pubdata/ipf_scrna-seq/ipf_obj_merged.qc_filtered.QConly_DFfiltered_singlets.rds"
OUT_DF_SUMMARY      <- "/data1/pubdata/ipf_scrna-seq/ipf_obj_merged.qc_filtered.QConly_DF_summary.tsv"
PLOT_OUT            <- "/data1/pubdata/ipf_scrna-seq/doublet_visualization.png"

PCS_USE    <- 1:20
RESOLUTION <- 0.8
EXP_RATE   <- 0.075
PN         <- 0.25


# =========================
# Helper: The "xtfrm" Killer
# =========================
fix_meta_atomic <- function(seu) {
  md <- as.data.frame(seu@meta.data)
  for (nm in colnames(md)) {
    if (is.list(md[[nm]]) || is.matrix(md[[nm]]) || is.data.frame(md[[nm]])) {
      md[[nm]] <- as.vector(unlist(md[[nm]]))
    }
  }
  seu@meta.data <- md
  return(seu)
}

# =========================
# DoubletFinder Core Function
# =========================
run_doubletfinder_sct <- function(o, npcs=30, pcs_use=1:20, resolution=0.8, pN=0.25, exp_rate=0.075) {
  # 1. SCTransform 메모리 절약 모드 활성화
  o <- SCTransform(o, verbose = FALSE, 
                   vst.flavor = "v2", 
                   conserve.memory = TRUE, # 중간 데이터 저장 억제
                   return.only.var.genes = TRUE) # 분석에 필요한 유전자만 유지
  
  o <- RunPCA(o, npcs = npcs, verbose = FALSE)
  
  # 2. pK Sweep 과정에서 가비지 컬렉션 강제 수행
  sweep.res <- paramSweep(o, PCs = pcs_use, sct = TRUE)
  gc() # 메모리 정리
  
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  best_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))  
  # nExp
  homotypic.prop <- modelHomotypic(o$seurat_clusters)
  nExp_adj <- round(round(exp_rate * ncol(o)) * (1 - homotypic.prop))
  
  # Run
  o <- doubletFinder(o, PCs = pcs_use, pN = pN, pK = best_pK, nExp = nExp_adj, sct = TRUE)
  
  class_col <- grep("^DF.classifications_", colnames(o@meta.data), value = TRUE) %>% tail(1)
  pann_col  <- grep("^pANN_", colnames(o@meta.data), value = TRUE) %>% tail(1)
  
  o$DF_class <- o@meta.data[[class_col]]
  o$DF_pANN  <- o@meta.data[[pann_col]]
  return(o)
}


# =========================
# MAIN - 병렬 실행 구간
# =========================
cat("[MAIN] Loading Object...\n")
obj <- readRDS(IN_RDS)
obj <- fix_meta_atomic(obj)

split_col <- if ("Sample" %in% colnames(obj@meta.data)) "Sample" else "orig.ident"
obj_list  <- SplitObject(obj, split.by = split_col)

cat(sprintf("[MAIN] Processing %d samples in parallel...\n", length(obj_list)))

# future_lapply를 사용하여 병렬 처리 수행
results_list <- future_lapply(names(obj_list), function(nm) {
  # 개별 샘플에 대해 DoubletFinder 실행
  # 내부 로그는 섞일 수 있으므로 temp_obj 생성에 집중
  t_obj <- run_doubletfinder_sct(obj_list[[nm]], pcs_use = 1:20, exp_rate = 0.075)
  t_obj <- fix_meta_atomic(t_obj)
  
  # 필요한 결과값만 리스트로 반환
  list(
    sample = nm,
    meta = t_obj@meta.data[, c("DF_class", "DF_pANN")],
    stats = data.frame(
      Sample = nm, 
      Total = ncol(t_obj),
      Singlets = as.integer(table(t_obj$DF_class)["Singlet"]),
      Doublets = as.integer(table(t_obj$DF_class)["Doublet"]),
      stringsAsFactors = FALSE
    )
  )
}, future.seed = TRUE)

# =========================
# 결과 통합
# =========================
cat("[MAIN] Merging results...\n")
final_meta_list <- lapply(results_list, function(x) x$meta)
summary_list <- lapply(results_list, function(x) x$stats)

final_meta_df <- do.call(rbind, final_meta_list)
obj$DF_class <- final_meta_df[rownames(obj@meta.data), "DF_class"]
obj$DF_pANN  <- final_meta_df[rownames(obj@meta.data), "DF_pANN"]

# # =========================
# # MAIN
# # =========================
# cat("[MAIN] Loading Object...\n")
# obj <- readRDS(IN_RDS)
# obj <- fix_meta_atomic(obj)
# 
# split_col <- if ("Sample" %in% colnames(obj@meta.data)) "Sample" else "orig.ident"
# obj_list  <- SplitObject(obj, split.by = split_col)
# 
# meta_to_add <- list()
# summary_list <- list()
# 
# for (nm in names(obj_list)) {
#   cat(sprintf("\n>>> Processing Sample: %s\n", nm))
#   temp_obj <- run_doubletfinder_sct(obj_list[[nm]], pcs_use = PCS_USE, resolution = RESOLUTION, exp_rate = EXP_RATE)
#   
#   temp_obj <- fix_meta_atomic(temp_obj)
#   meta_to_add[[nm]] <- temp_obj@meta.data[, c("DF_class", "DF_pANN")]
#   
#   stats <- table(temp_obj$DF_class)
#   summary_list[[nm]] <- data.frame(Sample=nm, Total=ncol(temp_obj), 
#                                    Singlets=as.integer(stats["Singlet"]), 
#                                    Doublets=as.integer(stats["Doublet"]), stringsAsFactors=F)
#   rm(temp_obj); gc()
# }
# 
# # Integration (Safe Mapping)
# final_meta_df <- do.call(rbind, meta_to_add)
# obj_barcodes <- rownames(obj@meta.data)
# final_meta_df <- final_meta_df[obj_barcodes, ]
# obj$DF_class <- final_meta_df$DF_class
# obj$DF_pANN  <- final_meta_df$DF_pANN

# =========================
# VISUALIZATION
# =========================
cat("\n[MAIN] Generating Visualizations...\n")
# Temporary processing for UMAP visualization
viz_obj <- obj
viz_obj <- NormalizeData(viz_obj, verbose = FALSE)
viz_obj <- FindVariableFeatures(viz_obj, verbose = FALSE)
viz_obj <- ScaleData(viz_obj, verbose = FALSE)
viz_obj <- RunPCA(viz_obj, verbose = FALSE)
viz_obj <- RunUMAP(viz_obj, dims = 1:20, verbose = FALSE)

# Plot 1: UMAP by Singlet/Doublet
p1 <- DimPlot(viz_obj, group.by = "DF_class", cols = c("Doublet" = "red", "Singlet" = "lightgrey"), pt.size = 0.5) + 
  ggtitle("DoubletFinder Classifications")

# Plot 2: Doublet Probability (pANN)
p2 <- FeaturePlot(viz_obj, features = "DF_pANN", pt.size = 0.5) + 
  scale_color_viridis_c(option = "plasma") + 
  ggtitle("Doublet Probability (pANN)")

ggsave(PLOT_OUT, p1 + p2, width = 14, height = 6)
cat(sprintf("[DONE] Visualization saved to: %s\n", PLOT_OUT))

# =========================
# SAVE DATA
# =========================
write.table(do.call(rbind, summary_list), OUT_DF_SUMMARY, sep="\t", quote=F, row.names=F)
saveRDS(obj, OUT_RDS_WITH_DF)
saveRDS(subset(obj, subset = DF_class == "Singlet"), OUT_RDS_DF_FILTERED)

cat("\n[COMPLETE] All tasks finished.\n")
# quick overall sanity
cat("\n[Overall DF_class]\n")
print(table(obj$DF_class, useNA = "ifany"))
cat("\n[Saved]\n- ", OUT_RDS_WITH_DF, "\n- ", OUT_DF_SUMMARY, "\n- ", OUT_RDS_DF_FILTERED, "\n")
