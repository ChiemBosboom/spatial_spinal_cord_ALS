#!/usr/bin/env Rscript

# 1. Setup Logging & Packages
log_file <- snakemake@log[[1]]
if (!is.null(log_file) && log_file != "") {
  log_con <- file(log_file, open = "wt")
  sink(log_con)
  sink(log_con, type = "message")
}

suppressPackageStartupMessages({
  library(TESSERA)
  library(fdrtool)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(Matrix)
  library(scales)
})

# Inputs, Outputs & Parameters
fits_rds      <- snakemake@input[["tessera_fits"]]
data_rds      <- snakemake@input[["tessera_data"]]
out_all_csv   <- snakemake@output[["results_all"]]
out_sig_csv   <- snakemake@output[["results_sig"]]
plots_dir     <- snakemake@output[["plots_dir"]]
contrasts_cfg <- snakemake@params[["contrasts"]]
fdr_thresh    <- as.numeric(snakemake@params[["fdr_threshold"]])

cat(sprintf("[%s] Starting TESSERA inference (FDR < %.2f | %d contrasts)...\n", 
            Sys.time(), fdr_thresh, length(contrasts_cfg)))

if (length(contrasts_cfg) == 0) {
  stop("FATAL: No contrasts provided in config.yaml.")
}

# -----------------------------------------------------------------------------
# 2. Load Models & Metadata
# -----------------------------------------------------------------------------
TESSERA_fits <- readRDS(fits_rds)
TESSERA_data <- readRDS(data_rds)

first_fit  <- TESSERA_fits[[1]]
coef_names <- names(first_fit$beta_hat)
cat("Model fixed effects:\n", paste("  -", coef_names, collapse = "\n"), "\n")

counts_list   <- TESSERA_data$counts_list
lib_size_list <- TESSERA_data$library_size_list
cov_list      <- TESSERA_data$covariates_list

# Extract sample info directly assuming "condition"
sample_info <- data.frame(
  sample    = names(counts_list),
  condition = vapply(names(counts_list), function(s) as.character(cov_list[[s]]$condition[1]), character(1)),
  stringsAsFactors = FALSE
) %>% arrange(condition, sample)

unique_conds <- unique(sample_info$condition)
cond_palette <- setNames(scales::hue_pal()(length(unique_conds)), unique_conds)

# -----------------------------------------------------------------------------
# 3. Dynamic Contrast Matrix Construction
# -----------------------------------------------------------------------------
parse_contrast_vector <- function(contrast_str, coef_names) {
  clean_coefs <- gsub(":", "__", coef_names)
  clean_str   <- gsub(":", "__", contrast_str)
  
  env <- new.env(parent = baseenv())
  p   <- length(coef_names)
  
  expr_parsed   <- parse(text = clean_str)
  contrast_vars <- all.vars(expr_parsed)
  
  for (cv in contrast_vars) {
    match_idx <- which(clean_coefs == cv)
    if (length(match_idx) == 0) match_idx <- which(tolower(clean_coefs) == tolower(cv))
    
    if (length(match_idx) == 1) {
      basis_vec <- numeric(p)
      basis_vec[match_idx] <- 1
      assign(cv, basis_vec, envir = env)
    } else {
      stop(sprintf("Variable '%s' in contrast '%s' could not be resolved.", cv, contrast_str))
    }
  }
  
  parsed_expr <- eval(expr_parsed, envir = env)
  names(parsed_expr) <- coef_names
  return(parsed_expr)
}

contrast_mat <- do.call(rbind, lapply(names(contrasts_cfg), function(c_name) {
  parse_contrast_vector(contrasts_cfg[[c_name]], coef_names)
}))
rownames(contrast_mat) <- names(contrasts_cfg)
colnames(contrast_mat) <- coef_names

# -----------------------------------------------------------------------------
# 4. Wald Tests & Per-Contrast Empirical Null Estimation
# -----------------------------------------------------------------------------
cat("Calculating Wald statistics across all genes...\n")
wald_list <- lapply(names(TESSERA_fits), function(g_name) {
  tryCatch({
    res <- TESSERA::calc_Wald_statistics(TESSERA_fits[[g_name]], contrast_mat)
    res$contrast_description <- rownames(res)
    res$gene <- g_name
    return(res)
  }, error = function(e) {
    warning(sprintf("Wald test failed for gene '%s': %s", g_name, conditionMessage(e)))
    return(NULL)
  })
})

wald_df <- do.call(rbind, wald_list)
wald_df <- wald_df[!is.na(wald_df$wald_stat_t) & is.finite(wald_df$wald_stat_t), ]
wald_df$wald_stat_chi2 <- wald_df$wald_stat_t^2

# Initialize inference tracking columns
wald_df$null_fit_status <- NA_character_
wald_df$null_scale      <- NA_real_
wald_df$null_shift      <- NA_real_
wald_df$null_threshold  <- NA_real_
wald_df$pval            <- NA_real_
wald_df$padj            <- NA_real_

cat("Estimating null distribution and computing p-values per contrast...\n")
for (cd in unique(wald_df$contrast_description)) {
  idx_cd     <- which(wald_df$contrast_description == cd)
  chi2_stats <- wald_df$wald_stat_chi2[idx_cd]
  t_stats    <- wald_df$wald_stat_t[idx_cd]
  
  # --- TIER 1: TESSERA select_Wald_threshold ---
  thresh_fit <- tryCatch({
    TESSERA::select_Wald_threshold(
      wald_stats       = chi2_stats,
      quantile_spacing = 0.05,  
      metric           = "Raw_MSE"
    )
  }, error = function(e) {
    cat(sprintf("  [%s] TESSERA thresholding failed (%s).\n", cd, conditionMessage(e)))
    return(NULL)
  })
  
  if (!is.null(thresh_fit) && !is.na(thresh_fit$chi2_params["scale"]) && thresh_fit$chi2_params["scale"] > 0) {
    scale_val <- as.numeric(thresh_fit$chi2_params["scale"])
    shift_val <- as.numeric(thresh_fit$chi2_params["shift"])
    cutoff    <- as.numeric(thresh_fit$threshold)
    
    cat(sprintf("  [%s] EMPIRICAL NULL (TESSERA): Scale=%.2f, Shift=%.4e, Cutoff=%.2f\n", 
                cd, scale_val, shift_val, cutoff))
    
    q_stat <- chi2_stats / scale_val
    wald_df$pval[idx_cd]            <- pchisq(q_stat, df = 1, ncp = shift_val, lower.tail = FALSE)
    wald_df$null_fit_status[idx_cd] <- "empirical_tessera"
    wald_df$null_scale[idx_cd]      <- scale_val
    wald_df$null_shift[idx_cd]      <- shift_val
    wald_df$null_threshold[idx_cd]  <- cutoff

  } else {
    # --- TIER 2: fdrtool FNDR Empirical Null  ---
    cat(sprintf("  [%s] Attempting empirical null calibration via fdrtool...\n", cd))
    
    fdr_res <- tryCatch({
      suppressMessages(
        fdrtool::fdrtool(
          x             = t_stats,
          statistic     = "normal",
          plot          = FALSE,
          cutoff.method = "fndr",
          verbose       = FALSE
        )
      )
    }, error = function(e) {
      cat(sprintf("  [%s] fdrtool failed (%s).\n", cd, conditionMessage(e)))
      return(NULL)
    })
    
    if (!is.null(fdr_res)) {
      sigma0 <- fdr_res$param[1, "sd"]
      cat(sprintf("  [%s] EMPIRICAL NULL (fdrtool): Null SD (sigma0)=%.2f (Equivalent Scale=%.2f)\n", 
                  cd, sigma0, sigma0^2))
      
      wald_df$pval[idx_cd]            <- fdr_res$pval
      wald_df$null_fit_status[idx_cd] <- "empirical_fdrtool"
      wald_df$null_scale[idx_cd]      <- sigma0^2
      wald_df$null_shift[idx_cd]      <- 0
      wald_df$null_threshold[idx_cd]  <- NA_real_

    } else {
      # --- TIER 3: Asymptotic Fallback (Last Resort) ---
      cat(sprintf("  [%s] WARNING: Both empirical methods failed. Using ASYMPTOTIC fallback.\n", cd))
      wald_df$pval[idx_cd]            <- pchisq(chi2_stats, df = 1, lower.tail = FALSE)
      wald_df$null_fit_status[idx_cd] <- "asymptotic_fallback"
    }
  }
  
  # Standard BH multiple testing correction
  wald_df$padj[idx_cd] <- p.adjust(wald_df$pval[idx_cd], method = "BH")
}

# Negative log10 padj capped at machine minimum to prevent infinite volcano coordinates
wald_df$neg_log10_padj <- -log10(pmax(wald_df$padj, 1e-300))
wald_df$is_significant <- !is.na(wald_df$padj) & (wald_df$padj < fdr_thresh)
wald_df$direction      <- ifelse(!wald_df$is_significant, "Not Significant",
                                 ifelse(wald_df$contrast_val > 0, "Up", "Down"))

# Asymptotic baseline reference columns
wald_df$pval_asymptotic <- pchisq(wald_df$wald_stat_chi2, df = 1, lower.tail = FALSE)
wald_df$padj_asymp      <- p.adjust(wald_df$pval_asymptotic, method = "BH")

# -----------------------------------------------------------------------------
# 5. Fast Gene Mean Expression & Scale Conversion
# -----------------------------------------------------------------------------
cat("Calculating gene-level average expression...\n")
unique_genes  <- unique(wald_df$gene)
all_lib_sizes <- unlist(lapply(lib_size_list, as.numeric))
total_spots   <- length(all_lib_sizes)

# Fast vectorized mean expression across spots without normalizing full 2000-gene matrices
total_counts <- setNames(numeric(length(unique_genes)), unique_genes)
for (s in names(counts_list)) {
  sub_counts <- counts_list[[s]][intersect(rownames(counts_list[[s]]), unique_genes), , drop = FALSE]
  total_counts[rownames(sub_counts)] <- total_counts[rownames(sub_counts)] + Matrix::rowSums(sub_counts)
}
mean_expr_df <- data.frame(
  gene            = names(total_counts),
  mean_expression = as.numeric(total_counts / max(total_spots, 1)),
  log2_mean_expr  = log2(as.numeric(total_counts / max(total_spots, 1)) + 1),
  stringsAsFactors = FALSE
)

wald_df <- dplyr::left_join(wald_df, mean_expr_df, by = "gene")
colnames(wald_df)[colnames(wald_df) == "contrast_val"] <- "estimate_ln"
colnames(wald_df)[colnames(wald_df) == "contrast_se"]  <- "std_error_ln"

# Convert GLMM natural log estimates to standard log2 fold change
wald_df$log2FC         <- wald_df$estimate_ln / log(2)
wald_df$std_error_log2 <- wald_df$std_error_ln / log(2)

# Top 5 genes for volcano annotation
wald_df$label_gene <- NA_character_
for (cd in unique(wald_df$contrast_description)) {
  idx_c    <- which(wald_df$contrast_description == cd)
  df_c     <- wald_df[idx_c, ]
  top_hits <- bind_rows(
    df_c %>% filter(is_significant & log2FC > 0) %>% arrange(padj, desc(abs(log2FC))) %>% head(5),
    df_c %>% filter(is_significant & log2FC < 0) %>% arrange(padj, desc(abs(log2FC))) %>% head(5)
  )
  flag_idx <- idx_c[wald_df$gene[idx_c] %in% top_hits$gene]
  wald_df$label_gene[flag_idx] <- wald_df$gene[flag_idx]
}

# -----------------------------------------------------------------------------
# 6. Export Results Tables
# -----------------------------------------------------------------------------
dir.create(dirname(out_all_csv), recursive = TRUE, showWarnings = FALSE)

out_cols <- c("gene", "contrast_description", "log2FC", "std_error_log2", 
              "estimate_ln", "std_error_ln", "wald_stat_t", "wald_stat_chi2", 
              "mean_expression", "log2_mean_expr", "null_fit_status", 
              "null_scale", "null_shift", "null_threshold", "pval", "padj", 
              "neg_log10_padj", "pval_asymptotic", "padj_asymp", "is_significant", "direction")

wald_out <- wald_df[, intersect(out_cols, colnames(wald_df))]
write.csv(wald_out, file = out_all_csv, row.names = FALSE)
cat("Full results table saved to:", out_all_csv, "\n")

sig_out <- wald_out %>% filter(is_significant == TRUE) %>% arrange(padj)
write.csv(sig_out, file = out_sig_csv, row.names = FALSE)
cat(sprintf("Significant results saved to: %s (%d hits at FDR < %.2f)\n", 
            out_sig_csv, nrow(sig_out), fdr_thresh))


# -----------------------------------------------------------------------------
# 7. Visualizations
# -----------------------------------------------------------------------------
volcano_dir <- file.path(plots_dir, "volcano")
ma_dir      <- file.path(plots_dir, "ma")
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ma_dir, recursive = TRUE, showWarnings = FALSE)

color_palette <- c("Up" = "#D73027", "Down" = "#4575B4", "Not Significant" = "#B0B0B0")

theme_boxed <- theme_bw(base_size = 12) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "#F2F2F2", color = "black", linewidth = 0.8),
    strip.text       = element_text(face = "bold", size = 11),
    axis.title       = element_text(face = "bold", size = 11),
    axis.text        = element_text(color = "black", size = 10),
    legend.position  = "top",
    legend.title     = element_text(face = "bold", size = 10),
    plot.title       = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle    = element_text(face = "italic", size = 10, hjust = 0.5, color = "grey30")
  )

# --- A. Moran's I Residual QC Boxplot ---
perf_list <- lapply(TESSERA_fits, function(x) x$performanceSummary)
perf_df   <- do.call(rbind, perf_list)

if (!is.null(perf_df) && all(c("Moran_counts", "Moran_residuals") %in% colnames(perf_df))) {
  perf_df <- perf_df %>% filter(!is.na(Moran_counts) & !is.na(Moran_residuals))
  moran_tidy <- data.frame(
    Sample = rep(perf_df$sample, 2),
    Type   = rep(c("Raw Counts", "Fitted Residuals"), each = nrow(perf_df)),
    MoranI = c(perf_df$Moran_counts, perf_df$Moran_residuals)
  )
  
  p_moran <- ggplot(moran_tidy, aes(x = Sample, y = MoranI, fill = Type)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.85, width = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_fill_manual(values = c("Raw Counts" = "#FDAE61", "Fitted Residuals" = "#2B83BA")) +
    theme_boxed +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Spatial Autocorrelation QC", x = "Sample", y = "Moran's I", fill = "Metric")
  
  ggsave(file.path(plots_dir, "moran_qc.png"), plot = p_moran, width = 8, height = 5, dpi = 300)
}

# --- B. Detect Mode: Continuous Gradient vs Categorical Sub-Niche vs Pooled Niche ---
sample_covs1 <- cov_list[[1]]
cov_names    <- colnames(sample_covs1)

# Identify non-condition variables present in the model formula
target_covs <- cov_names[cov_names != "condition" & vapply(cov_names, function(v) any(grepl(v, coef_names)), logical(1))]

if (length(target_covs) == 0) {
  # POOLED NICHE MODE: Formula is '~ 0 + condition' (evaluating all filtered spots together)
  is_continuous_analysis <- FALSE
  is_pooled_niche        <- TRUE
  target_cov             <- "condition"
  heat_dir               <- file.path(plots_dir, "heatmap")
  dir.create(heat_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Analysis mode: [Pooled Niche Comparison] across filtered spots.\n")
} else {
  is_pooled_niche        <- FALSE
  target_cov             <- target_covs[1]
  is_continuous_analysis <- is.numeric(sample_covs1[[target_cov]])
  
  if (is_continuous_analysis) {
    gradient_dir <- file.path(plots_dir, "gradient_profiles")
    dir.create(gradient_dir, recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("Analysis mode: [Continuous Gradient] on '%s'.\n", target_cov))
  } else {
    heat_dir <- file.path(plots_dir, "heatmap")
    dir.create(heat_dir, recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("Analysis mode: [Categorical Sub-Niche] on '%s'.\n", target_cov))
  }
}

# --- Helper: Normalize only the requested subset of genes on-the-fly ---
get_normalized_submat <- function(samples, genes) {
  mean_lib <- mean(all_lib_sizes, na.rm = TRUE)
  sub_list <- list()
  for (s in samples) {
    mat  <- counts_list[[s]]
    g_in <- intersect(rownames(mat), genes)
    lib  <- as.numeric(lib_size_list[[s]])
    lib[lib == 0 | is.na(lib)] <- 1
    norm_mat <- mat[g_in, , drop = FALSE] %*% Matrix::Diagonal(x = mean_lib / lib)
    sub_list[[s]] <- norm_mat
  }
  return(sub_list)
}

# --- C. Dedicated Plots per Contrast ---
for (cd in unique(wald_df$contrast_description)) {
  clean_cd <- gsub("[^A-Za-z0-9_.-]", "_", cd)
  df_c     <- wald_df[wald_df$contrast_description == cd, ]
  
  sub_label <- if (df_c$null_fit_status[1] == "asymptotic_fallback") {
    "Note: Empirical null failed; asymptotic p-values shown"
  } else NULL

  # 1. Dedicated Volcano Plot
  p_volc <- ggplot(df_c, aes(x = log2FC, y = neg_log10_padj)) +
    geom_point(data = filter(df_c, !is_significant),
               color = color_palette["Not Significant"], alpha = 0.35, size = 1.3)
  
  # Only plot significant points layer if hits exist to avoid ggplot scale warnings
  if (any(df_c$is_significant)) {
    p_volc <- p_volc +
      geom_point(data = filter(df_c, is_significant),
                 aes(color = direction), alpha = 0.85, size = 1.8)
  }
  
  if (any(!is.na(df_c$label_gene))) {
    p_volc <- p_volc +
      geom_point(data = filter(df_c, !is.na(label_gene)),
                 shape = 21, color = "black", fill = NA, size = 2.4, stroke = 0.8) +
      geom_text_repel(aes(label = label_gene),
                      data = filter(df_c, !is.na(label_gene)),
                      size = 3.6, fontface = "bold",
                      box.padding = 0.5, point.padding = 0.3,
                      max.overlaps = Inf, min.segment.length = 0,
                      segment.color = "grey40", segment.size = 0.3,
                      show.legend = FALSE)
  }

  p_volc <- p_volc +
    geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_color_manual(values = color_palette, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    theme_boxed +
    labs(
      title    = cd,
      subtitle = sub_label,
      x        = expression(bold(log[2]~"Fold Change")), 
      y        = expression(bold(-log[10]~"(adj. P-value)")), 
      color    = "Direction"
    )
  ggsave(file.path(volcano_dir, paste0("volcano_", clean_cd, ".png")), plot = p_volc, width = 7, height = 6, dpi = 300)

  # 2. Dedicated MA Plot
  p_ma <- ggplot(df_c, aes(x = log2_mean_expr, y = log2FC)) +
    geom_point(data = filter(df_c, !is_significant),
               color = color_palette["Not Significant"], alpha = 0.35, size = 1.3)
  
  if (any(df_c$is_significant)) {
    p_ma <- p_ma +
      geom_point(data = filter(df_c, is_significant),
                 aes(color = direction), alpha = 0.85, size = 1.8)
  }

  p_ma <- p_ma +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_color_manual(values = color_palette, drop = FALSE) +
    theme_boxed +
    labs(
      title = cd, 
      x     = expression(bold("Average Expression ("*log[2]*")")), 
      y     = expression(bold(log[2]~"Fold Change")), 
      color = "Direction"
    )
  ggsave(file.path(ma_dir, paste0("ma_", clean_cd, ".png")), plot = p_ma, width = 7, height = 5.5, dpi = 300)

  # 3. Third Plot: Gradient Profile or Categorical Heatmap
  active_coefs <- coef_names[contrast_mat[cd, ] != 0]
  
  # Determine active conditions
  active_conds <- unique(sample_info$condition[sapply(sample_info$condition, function(cond) {
    any(grepl(paste0("condition", cond), active_coefs, ignore.case = TRUE))
  })])
  if (length(active_conds) == 0) active_conds <- unique(sample_info$condition)
  active_samples <- sample_info %>% filter(condition %in% active_conds)

  if (is_continuous_analysis) {
    # -------------------------------------------------------------------------
    # MODE A: 10-Bin Mean +/- SE Gradient Profiles (Top 4 Genes)
    # -------------------------------------------------------------------------
    top4_genes <- df_c %>% arrange(padj, desc(abs(log2FC))) %>% head(4) %>% pull(gene)
    sub_norm   <- get_normalized_submat(active_samples$sample, top4_genes)
    
    spot_data_list <- list()
    for (s in active_samples$sample) {
      s_cond   <- active_samples$condition[active_samples$sample == s]
      grad_val <- as.numeric(cov_list[[s]][[target_cov]])
      norm_mat <- sub_norm[[s]][top4_genes, , drop = FALSE]
      log_expr <- as.matrix(log2(norm_mat + 1))
      
      for (g in top4_genes) {
        spot_data_list[[paste0(s, "_", g)]] <- data.frame(
          gene       = g,
          sample     = s,
          condition  = s_cond,
          gradient   = grad_val,
          expression = log_expr[g, ],
          stringsAsFactors = FALSE
        )
      }
    }
    
    df_spots  <- do.call(rbind, spot_data_list)
    g_range   <- range(df_spots$gradient, na.rm = TRUE)
    breaks    <- seq(g_range[1], g_range[2], length.out = 11)
    bin_width <- breaks[2] - breaks[1]
    midpoints <- (breaks[-1] + breaks[-length(breaks)]) / 2
    
    df_spots$bin_idx <- cut(df_spots$gradient, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    df_spots$bin_x   <- midpoints[df_spots$bin_idx]
    
    bin_summary <- df_spots %>%
      filter(!is.na(bin_x)) %>%
      group_by(gene, condition, bin_x) %>%
      summarise(
        mean_expr = mean(expression, na.rm = TRUE),
        se_expr   = ifelse(n() > 1, sd(expression, na.rm = TRUE) / sqrt(n()), 0),
        n_spots   = n(),
        .groups   = "drop"
      ) %>%
      mutate(
        ymin = mean_expr - se_expr,
        ymax = mean_expr + se_expr,
        gene = factor(gene, levels = top4_genes),
        condition = factor(condition, levels = active_conds)
      )
    
    pd <- if (length(active_conds) > 1) position_dodge(width = bin_width * 0.25) else position_identity()
    dynamic_x_label <- tools::toTitleCase(gsub("_", " ", target_cov))
    
    p_grad <- ggplot(bin_summary, aes(x = bin_x, y = mean_expr, color = condition, group = condition)) +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = bin_width * 0.3, position = pd, alpha = 0.75, linewidth = 0.6) +
      geom_line(position = pd, linewidth = 1.0) +
      geom_point(position = pd, size = 2.4) +
      scale_color_manual(values = cond_palette) +
      facet_wrap(~gene, scales = "free_y", ncol = 2) +
      theme_boxed +
      labs(
        title = cd,
        x     = dynamic_x_label,
        y     = expression(bold("Normalized Expression ("*log[2]*")")),
        color = "Condition"
      )
    
    grad_file <- file.path(gradient_dir, paste0("gradient_profile_", clean_cd, ".png"))
    ggsave(grad_file, plot = p_grad, width = 9.5, height = 7.5, dpi = 300)

  } else {
    # -------------------------------------------------------------------------
    # MODE B: Niche Heatmap 
    # -------------------------------------------------------------------------
    top_sig <- df_c %>% arrange(padj, desc(abs(log2FC))) %>% head(30)
    if (nrow(top_sig) == 0) next
    
    top_sig     <- top_sig %>% arrange(log2FC)
    gene_levels <- top_sig$gene
    sub_norm    <- get_normalized_submat(active_samples$sample, gene_levels)
    
    cell_list <- list()

    if (is_pooled_niche) {
      # 1. POOLED NICHE: Average across all filtered spots per sample
      cat(sprintf("    Generating pooled niche heatmap for [%s]...\n", 
                  paste(active_conds, collapse = ", ")))
      
      for (s in active_samples$sample) {
        s_cond   <- active_samples$condition[active_samples$sample == s]
        norm_mat <- sub_norm[[s]][gene_levels, , drop = FALSE]
        mean_exp <- Matrix::rowMeans(norm_mat)
        
        cell_list[[s]] <- data.frame(
          gene       = gene_levels,
          sample     = s,
          condition  = s_cond,
          level      = "Microenvironment",
          expression = log2(mean_exp + 1),
          stringsAsFactors = FALSE
        )
      }
    } else {
      # 2. SUB-NICHE: Split by categorical levels (e.g., True/False)
      all_levels <- unique(as.character(sample_covs1[[target_cov]]))
      active_levels <- all_levels[sapply(all_levels, function(lev) {
        any(grepl(paste0(target_cov, lev), active_coefs, ignore.case = TRUE))
      })]
      if (length(active_levels) == 0) active_levels <- all_levels
      
      cat(sprintf("    Generating sub-niche heatmap for [%s] across levels [%s]...\n", 
                  paste(active_conds, collapse = ", "), paste(active_levels, collapse = ", ")))
      
      for (s in active_samples$sample) {
        s_cond   <- active_samples$condition[active_samples$sample == s]
        cov_s    <- as.character(cov_list[[s]][[target_cov]])
        norm_mat <- sub_norm[[s]][gene_levels, , drop = FALSE]
        
        for (lev in active_levels) {
          spots_idx <- which(tolower(cov_s) == tolower(lev))
          if (length(spots_idx) == 0) next
          
          mean_lev <- Matrix::rowMeans(norm_mat[, spots_idx, drop = FALSE])
          
          cell_list[[paste0(s, "_", lev)]] <- data.frame(
            gene       = gene_levels,
            sample     = s,
            condition  = s_cond,
            level      = lev,
            expression = log2(mean_lev + 1),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    df_raw <- do.call(rbind, cell_list)
    
    # Compute Balanced Z-scores
    df_heat <- df_raw %>%
      group_by(gene, level) %>%
      mutate(
        balanced_center = mean(tapply(expression, condition, mean)),
        sd_val          = sd(expression),
        score           = if (!is.na(sd_val[1]) && sd_val[1] > 0) {
          (expression - balanced_center) / sd_val[1]
        } else {
          0
        }
      ) %>%
      ungroup()
    
    # Order factors
    if (!is_pooled_niche) {
      lev_unique <- unique(df_heat$level)
      if (all(c("true", "false") %in% tolower(lev_unique))) {
        lev_order <- lev_unique[order(tolower(lev_unique) != "true")]
      } else {
        lev_order <- sort(lev_unique)
      }
      df_heat$level <- factor(df_heat$level, levels = lev_order)
    }

    df_heat$condition <- factor(df_heat$condition, levels = active_conds)
    df_heat$gene      <- factor(df_heat$gene, levels = gene_levels)
    df_heat$sample    <- factor(df_heat$sample, levels = rev(unique(df_heat$sample)))
    
    # Determine faceting formula
    has_multiple_levels <- (!is_pooled_niche) && (length(unique(df_heat$level)) > 1)
    facet_formula       <- if (has_multiple_levels) level + condition ~ . else condition ~ .
    
    subtitle_text <- if (is_pooled_niche) {
      ""
    } else if (!has_multiple_levels) {
      paste0(sub_label, if (!is.null(sub_label)) " | " else "", target_cov, ": ", unique(df_heat$level))
    } else {
      sub_label
    }

    p_heat <- ggplot(df_heat, aes(x = gene, y = sample, fill = score)) +
      geom_tile(color = "white", linewidth = 0.5) +
      facet_grid(facet_formula, scales = "free_y", space = "free_y") +
      scale_fill_gradient2(
        low      = "#4575B4", 
        mid      = "white", 
        high     = "#D73027", 
        midpoint = 0, 
        limits   = c(-2.5, 2.5),
        oob      = scales::squish,
        name     = "Balanced Z-Score",
        guide    = guide_colorbar(
          barwidth       = unit(5, "cm"),
          barheight      = unit(0.35, "cm"),
          title.position = "top",
          title.hjust    = 0.5
        )
      ) +
      theme_boxed +
      theme(
        axis.text.x      = element_text(angle = 45, hjust = 1, face = "bold.italic", size = 9),
        axis.text.y      = element_text(face = "bold", size = 9, color = "black"),
        strip.background = element_rect(fill = "#EFEFEF", color = "black", linewidth = 0.6),
        strip.text.y     = element_text(face = "bold", size = 9, angle = 0),
        panel.spacing.y  = unit(0.3, "lines"),
        legend.text      = element_text(size = 8),
        legend.title     = element_text(face = "bold", size = 9, hjust = 0.5)
      ) +
      labs(
        title    = cd, 
        subtitle = subtitle_text,
        x        = NULL, 
        y        = NULL
      )
    
    n_sample_rows <- nrow(distinct(df_heat, level, condition, sample))
    plot_height   <- max(3.0, 0.35 * n_sample_rows + 2.0)
    
    heat_file <- file.path(heat_dir, paste0("heatmap_", clean_cd, ".png"))
    ggsave(heat_file, plot = p_heat, width = 10.5, height = plot_height, dpi = 300)
    cat("    Heatmap saved to:", heat_file, "\n")
  }
}

cat(sprintf("[%s] TESSERA inference and plotting completed successfully!\n", Sys.time()))