#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# infer_tessera.R
# Conducts linear contrasts, empirical null distribution fitting, FDR correction,
# and generates dedicated publication-ready figures for every comparison:
#   - Volcano Plot (per contrast)
#   - MA Plot (per contrast)
#   - Moran's I Residual QC Plot (overall)
#   - Categorical Analysis: Sample x Compartment Expression Heatmap
#   - Continuous Analysis: 10-Bin Mean +/- SE Gradient Profiles (Top 4 Genes, Active Conditions Only)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 1. Setup Logging & Packages
# -----------------------------------------------------------------------------
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
})

# Inputs & Outputs
fits_rds    <- snakemake@input[["tessera_fits"]]
data_rds    <- snakemake@input[["tessera_data"]]

out_all_csv <- snakemake@output[["results_all"]]
out_sig_csv <- snakemake@output[["results_sig"]]
plots_dir   <- snakemake@output[["plots_dir"]]

# Parameters
contrasts_cfg <- snakemake@params[["contrasts"]]
fdr_thresh    <- as.numeric(snakemake@params[["fdr_threshold"]])

cat(sprintf("[%s] Starting TESSERA statistical inference...\n", Sys.time()))
cat(sprintf("FDR Threshold: %.2f | Number of configured contrasts: %d\n", 
            fdr_thresh, length(contrasts_cfg)))

if (length(contrasts_cfg) == 0) {
  stop("FATAL: No contrasts were provided in config.yaml for this analysis.")
}

# -----------------------------------------------------------------------------
# 2. Load Fitted Models and Data Container
# -----------------------------------------------------------------------------
cat("Loading fitted models from:", fits_rds, "\n")
TESSERA_fits <- readRDS(fits_rds)
cat(sprintf("Loaded fits for %d genes.\n", length(TESSERA_fits)))

cat("Loading data container from:", data_rds, "\n")
TESSERA_data <- readRDS(data_rds)

first_fit  <- TESSERA_fits[[1]]
coef_names <- names(first_fit$beta_hat)
cat("Model fixed effect parameters:\n")
cat(paste("  -", coef_names, collapse = "\n"), "\n")

# -----------------------------------------------------------------------------
# 3. Dynamic Contrast Matrix Construction
# -----------------------------------------------------------------------------
parse_contrast_vector <- function(contrast_str, coef_names) {
  clean_coefs <- gsub(":", "__", coef_names)
  clean_str   <- gsub(":", "__", contrast_str)
  
  env <- new.env(parent = baseenv())
  p <- length(coef_names)
  
  expr_parsed <- tryCatch({
    parse(text = clean_str)
  }, error = function(e) {
    stop(sprintf("Syntax error in contrast: '%s'. Error: %s", contrast_str, conditionMessage(e)))
  })
  
  contrast_vars <- all.vars(expr_parsed)
  
  for (cv in contrast_vars) {
    match_idx <- which(clean_coefs == cv)
    if (length(match_idx) == 0) {
      match_idx <- which(tolower(clean_coefs) == tolower(cv))
    }
    
    if (length(match_idx) == 1) {
      basis_vec <- numeric(p)
      basis_vec[match_idx] <- 1
      assign(cv, basis_vec, envir = env)
    } else if (length(match_idx) == 0) {
      stop(sprintf("Variable '%s' in contrast '%s' does not match any coefficient.\nAvailable:\n  %s",
                   cv, contrast_str, paste(coef_names, collapse = "\n  ")))
    } else {
      stop(sprintf("Variable '%s' in contrast '%s' is ambiguous.", cv, contrast_str))
    }
  }
  
  parsed_expr <- eval(expr_parsed, envir = env)
  if (length(parsed_expr) != p || !is.numeric(parsed_expr)) {
    stop(sprintf("Contrast '%s' did not evaluate to a valid numeric vector.", contrast_str))
  }
  
  names(parsed_expr) <- coef_names
  return(parsed_expr)
}

cat("Constructing contrast matrix...\n")
contrast_vec_list <- lapply(names(contrasts_cfg), function(c_name) {
  parse_contrast_vector(contrasts_cfg[[c_name]], coef_names)
})
contrast_mat <- do.call(rbind, contrast_vec_list)
rownames(contrast_mat) <- names(contrasts_cfg)
colnames(contrast_mat) <- coef_names

cat("Contrast matrix preview:\n")
print(contrast_mat)

# -----------------------------------------------------------------------------
# 4. Compute Wald Statistics
# -----------------------------------------------------------------------------
cat("Calculating Wald statistics across all genes and contrasts...\n")
wald_list <- lapply(names(TESSERA_fits), function(g_name) {
  tryCatch({
    fit <- TESSERA_fits[[g_name]]
    res <- TESSERA::calc_Wald_statistics(fit, contrast_mat)
    res$contrast_description <- rownames(res)
    res$gene <- g_name
    return(res)
  }, error = function(e) {
    warning(sprintf("Wald calculation failed for gene '%s': %s", g_name, conditionMessage(e)))
    return(NULL)
  })
})

wald_df <- do.call(rbind, wald_list)
wald_df <- wald_df[!is.na(wald_df$wald_stat_t) & is.finite(wald_df$wald_stat_t), ]
cat(sprintf("Calculated Wald tests for %d gene-contrast comparisons.\n", nrow(wald_df)))

# -----------------------------------------------------------------------------
# 5. Empirical Null Distribution Estimation & Exact Log P-Values
# -----------------------------------------------------------------------------
cat("Fitting empirical null distributions...\n")
wald_df$wald_stat_chi2 <- wald_df$wald_stat_t^2

null_scale <- 1
null_shift <- 0
thresh_success <- FALSE

tryCatch({
  thresh_fit <- TESSERA::select_Wald_threshold(
    wald_stats       = wald_df$wald_stat_chi2,
    quantile_spacing = 0.01,
    metric           = "Raw_MSE"
  )
  null_scale     <- thresh_fit$chi2_params["scale"]
  null_shift     <- thresh_fit$chi2_params["shift"]
  thresh_success <- TRUE
  cat(sprintf("  TESSERA empirical null: Scale=%.4e, Shift=%.4e, Cutoff=%.2f\n", 
              null_scale, null_shift, thresh_fit$threshold))
}, error = function(e) {
  warning("select_Wald_threshold failed; falling back to asymptotic chi^2: ", conditionMessage(e))
})

# Exact log(p-value) prevents numerical underflow to 0
q_stat <- wald_df$wald_stat_chi2 / null_scale
wald_df$log_pval_tessera <- pchisq(q_stat, df = 1, ncp = null_shift, lower.tail = FALSE, log.p = TRUE)
wald_df$pval_tessera     <- exp(wald_df$log_pval_tessera)

# -----------------------------------------------------------------------------
# 6. Multiple Testing Correction in Log-Space
# -----------------------------------------------------------------------------
calc_log10_padj <- function(log_p) {
  m <- length(log_p)
  ord <- order(log_p)
  log_padj_sorted <- log_p[ord] + log(m) - log(seq_len(m))
  log_padj_sorted <- rev(cummin(rev(log_padj_sorted)))
  log_padj_sorted <- pmin(log_padj_sorted, 0)
  
  log_padj <- numeric(m)
  log_padj[ord] <- log_padj_sorted
  neg_log10_padj <- -log_padj / log(10)
  return(list(padj = exp(log_padj), neg_log10_padj = neg_log10_padj))
}

wald_df$padj_tessera   <- NA_real_
wald_df$neg_log10_padj <- NA_real_

for (cd in unique(wald_df$contrast_description)) {
  idx <- which(wald_df$contrast_description == cd)
  adj_res <- calc_log10_padj(wald_df$log_pval_tessera[idx])
  wald_df$padj_tessera[idx]   <- adj_res$padj
  wald_df$neg_log10_padj[idx] <- adj_res$neg_log10_padj
}

wald_df$is_significant <- !is.na(wald_df$padj_tessera) & (wald_df$padj_tessera < fdr_thresh)
wald_df$direction      <- ifelse(!wald_df$is_significant, "Not Significant",
                                 ifelse(wald_df$contrast_val > 0, "Up", "Down"))

# Asymptotic p-values for reference
wald_df$pval_asymptotic <- pchisq(wald_df$wald_stat_chi2, df = 1, lower.tail = FALSE)
wald_df$padj_asymp      <- p.adjust(wald_df$pval_asymptotic, method = "BH")

# -----------------------------------------------------------------------------
# 7. Extract Library-Size Normalized Expression per Spot
# -----------------------------------------------------------------------------
cat("Calculating library-size normalized expression per gene...\n")
unique_genes  <- unique(wald_df$gene)
counts_list   <- TESSERA_data$counts_list
lib_size_list <- TESSERA_data$library_size_list

all_lib_sizes <- unlist(lapply(lib_size_list, as.numeric))
mean_lib_size <- mean(all_lib_sizes, na.rm = TRUE)
total_spots   <- length(all_lib_sizes)
genes_in_rows <- any(unique_genes %in% rownames(counts_list[[1]]))

total_norm_counts    <- setNames(numeric(length(unique_genes)), unique_genes)
sample_norm_mat_list <- list()

for (s in names(counts_list)) {
  mat <- counts_list[[s]]
  lib <- as.numeric(lib_size_list[[s]])
  lib[lib == 0 | is.na(lib)] <- 1
  norm_factors <- mean_lib_size / lib
  
  if (genes_in_rows) {
    common_g <- intersect(rownames(mat), unique_genes)
    sub_mat  <- mat[common_g, , drop = FALSE]
    norm_mat <- sub_mat %*% Matrix::Diagonal(x = norm_factors)
    total_norm_counts[common_g] <- total_norm_counts[common_g] + Matrix::rowSums(norm_mat)
    sample_norm_mat_list[[s]]   <- norm_mat
  } else {
    common_g <- intersect(colnames(mat), unique_genes)
    sub_mat  <- mat[, common_g, drop = FALSE]
    norm_mat <- Matrix::Diagonal(x = norm_factors) %*% sub_mat
    total_norm_counts[common_g] <- total_norm_counts[common_g] + Matrix::colSums(norm_mat)
    sample_norm_mat_list[[s]]   <- t(norm_mat)
  }
}

mean_expr <- total_norm_counts / max(total_spots, 1)
mean_expr_df <- data.frame(
  gene            = names(mean_expr),
  mean_expression = as.numeric(mean_expr),
  log2_mean_expr  = log2(as.numeric(mean_expr) + 1),
  stringsAsFactors = FALSE
)

wald_df <- dplyr::left_join(wald_df, mean_expr_df, by = "gene")
colnames(wald_df)[colnames(wald_df) == "contrast_val"] <- "estimate"
colnames(wald_df)[colnames(wald_df) == "contrast_se"]  <- "std_error"

# -----------------------------------------------------------------------------
# 8. Identify Top 5 Up and Down Genes per Contrast (for Volcano Plot)
# -----------------------------------------------------------------------------
cat("Selecting top 5 up and down genes per contrast for labeling...\n")
wald_df$label_gene <- NA_character_

for (cd in unique(wald_df$contrast_description)) {
  idx_c <- which(wald_df$contrast_description == cd)
  df_c  <- wald_df[idx_c, ]
  
  top_up <- df_c %>%
    filter(is_significant & estimate > 0) %>%
    arrange(desc(neg_log10_padj), desc(estimate)) %>%
    head(5)
    
  top_down <- df_c %>%
    filter(is_significant & estimate < 0) %>%
    arrange(desc(neg_log10_padj), estimate) %>%
    head(5)
    
  labeled_genes_cd <- c(top_up$gene, top_down$gene)
  flag_idx <- idx_c[wald_df$gene[idx_c] %in% labeled_genes_cd]
  wald_df$label_gene[flag_idx] <- wald_df$gene[flag_idx]
}

# -----------------------------------------------------------------------------
# 9. Export Results Tables
# -----------------------------------------------------------------------------
dir.create(dirname(out_all_csv), recursive = TRUE, showWarnings = FALSE)

out_cols <- c("gene", "contrast_description", "estimate", "std_error", 
              "wald_stat_t", "wald_stat_chi2", "mean_expression", "log2_mean_expr",
              "pval_tessera", "padj_tessera", "neg_log10_padj", 
              "pval_asymptotic", "padj_asymp", "is_significant", "direction")

wald_out <- wald_df[, intersect(out_cols, colnames(wald_df))]
write.csv(wald_out, file = out_all_csv, row.names = FALSE)
cat("Full results table saved to:", out_all_csv, "\n")

sig_out <- wald_out %>% filter(is_significant == TRUE) %>% arrange(padj_tessera)
write.csv(sig_out, file = out_sig_csv, row.names = FALSE)
cat(sprintf("Significant results table saved to: %s (%d hits at FDR < %.2f)\n", 
            out_sig_csv, nrow(sig_out), fdr_thresh))

# -----------------------------------------------------------------------------
# 10. Publication-Ready Visualizations
# -----------------------------------------------------------------------------
volcano_dir <- file.path(plots_dir, "volcano")
ma_dir      <- file.path(plots_dir, "ma")

dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ma_dir, recursive = TRUE, showWarnings = FALSE)

color_palette <- c("Up" = "#D73027", "Down" = "#4575B4", "Not Significant" = "#B0B0B0")
cond_palette  <- c("ALS" = "#D73027", "CTRL" = "#4575B4")

theme_boxed <- theme_bw(base_size = 12) +
  theme(
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = "#F2F2F2", color = "black", linewidth = 0.8),
    strip.text         = element_text(face = "bold", size = 11),
    axis.title         = element_text(face = "bold", size = 11),
    axis.text          = element_text(color = "black", size = 10),
    legend.position    = "top",
    legend.title       = element_text(face = "bold", size = 10),
    plot.title         = element_text(face = "bold", size = 13, hjust = 0.5)
  )

# --- A. Moran's I QC Boxplot ---
cat("Generating Moran's I QC plot...\n")
perf_list <- lapply(TESSERA_fits, function(x) x$performanceSummary)
perf_df   <- do.call(rbind, perf_list)

if (!is.null(perf_df) && all(c("Moran_counts", "Moran_residuals") %in% colnames(perf_df))) {
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
    labs(title = "Spatial Autocorrelation (Moran's I)", x = "Sample", y = "Moran's I", fill = "Metric")
  
  ggsave(file.path(plots_dir, "moran_qc.png"), plot = p_moran, width = 8, height = 5, dpi = 300)
}

# --- Sample Info ---
sample_conditions <- sapply(names(counts_list), function(s) {
  if (!is.null(TESSERA_data$covariates_list[[s]]$condition)) {
    as.character(TESSERA_data$covariates_list[[s]]$condition[1])
  } else "Sample"
})

sample_info <- data.frame(
  sample    = names(counts_list),
  condition = sample_conditions,
  stringsAsFactors = FALSE
) %>% arrange(condition, sample)

# --- B. Detect Analysis Mode (Formula Variables Only) ---
sample_covs    <- TESSERA_data$covariates_list[[1]]
cov_candidates <- setdiff(colnames(sample_covs), c("condition", "sample_id", "barcode", "nested_id", "nUMI", "nGene"))

# Filter to variables that ACTUALLY appear in this model's coefficient names
formula_covs <- cov_candidates[sapply(cov_candidates, function(v) any(grepl(v, coef_names)))]
cat(sprintf("Variables detected in model formula: %s\n", paste(formula_covs, collapse = ", ")))

# Continuous analysis only if a numeric variable is part of the formula
is_continuous_analysis <- any(sapply(formula_covs, function(v) is.numeric(sample_covs[[v]])))

if (is_continuous_analysis) {
  continuous_cov_name <- formula_covs[which(sapply(formula_covs, function(v) is.numeric(sample_covs[[v]])))[1]]
  gradient_dir        <- file.path(plots_dir, "gradient_profiles")
  dir.create(gradient_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("Analysis mode: [Continuous Gradient] on '%s'. Generating 10-bin gradient profile plots.\n", continuous_cov_name))
} else {
  target_cov <- formula_covs[1]
  heat_dir   <- file.path(plots_dir, "heatmap")
  dir.create(heat_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("Analysis mode: [Categorical Niche] on '%s'. Generating Sample x Compartment heatmaps.\n", target_cov))
}

# --- C. Dedicated Plots per Comparison ---
cat("Generating dedicated figures per comparison...\n")

for (cd in unique(wald_df$contrast_description)) {
  clean_cd <- gsub("[^A-Za-z0-9_.-]", "_", cd)
  df_c     <- wald_df[wald_df$contrast_description == cd, ]
  
  cat(sprintf("  Processing contrast: %s\n", cd))

  # 1. Dedicated Volcano Plot
  p_volc <- ggplot(df_c, aes(x = estimate, y = neg_log10_padj)) +
    geom_point(data = filter(df_c, !is_significant),
               color = color_palette["Not Significant"], alpha = 0.35, size = 1.3) +
    geom_point(data = filter(df_c, is_significant),
               aes(color = direction), alpha = 0.85, size = 1.8) +
    geom_point(data = filter(df_c, !is.na(label_gene)),
               shape = 21, color = "black", fill = NA, size = 2.4, stroke = 0.8) +
    geom_text_repel(aes(label = label_gene),
                    data = filter(df_c, !is.na(label_gene)),
                    size = 3.6, fontface = "bold",
                    box.padding = 0.5, point.padding = 0.3,
                    max.overlaps = Inf, min.segment.length = 0,
                    segment.color = "grey40", segment.size = 0.3,
                    show.legend = FALSE) +
    geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_color_manual(values = color_palette) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    theme_boxed +
    labs(title = cd, x = "LogFC", y = "-log10(adj. P-value)", color = "Direction")
  
  ggsave(file.path(volcano_dir, paste0("volcano_", clean_cd, ".png")), plot = p_volc, width = 7, height = 6, dpi = 300)

  # 2. Dedicated MA Plot
  p_ma <- ggplot(df_c, aes(x = log2_mean_expr, y = estimate)) +
    geom_point(data = filter(df_c, !is_significant),
               color = color_palette["Not Significant"], alpha = 0.35, size = 1.3) +
    geom_point(data = filter(df_c, is_significant),
               aes(color = direction), alpha = 0.85, size = 1.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_color_manual(values = color_palette) +
    theme_boxed +
    labs(title = cd, x = "Average Expression (log2)", y = "LogFC", color = "Direction")
  
  ggsave(file.path(ma_dir, paste0("ma_", clean_cd, ".png")), plot = p_ma, width = 7, height = 5.5, dpi = 300)

  # 3. Third Plot: Continuous Gradient Profiles OR Categorical Heatmap
  active_coefs <- coef_names[contrast_mat[cd, ] != 0]
  
  # Determine which condition(s) are active in this contrast
  active_conds <- unique(sample_info$condition[sapply(sample_info$condition, function(cond) {
    any(grepl(paste0("condition", cond), active_coefs, ignore.case = TRUE))
  })])
  active_samples <- sample_info %>% filter(condition %in% active_conds)

  if (is_continuous_analysis) {
    # -------------------------------------------------------------------------
    # MODE A: 10-Bin Mean +/- SE Gradient Profiles (Top 4 Genes, Active Lines Only)
    # -------------------------------------------------------------------------
    top4_genes <- df_c %>% arrange(padj_tessera, desc(abs(estimate))) %>% head(4) %>% pull(gene)
    cat(sprintf("    Generating gradient profile for conditions [%s] (Top 4 genes: %s)...\n", 
                paste(active_conds, collapse = ", "), paste(top4_genes, collapse = ", ")))
    
    # Gather spot data only for the active samples/conditions
    spot_data_list <- list()
    for (s in active_samples$sample) {
      s_cond   <- active_samples$condition[active_samples$sample == s]
      grad_val <- as.numeric(TESSERA_data$covariates_list[[s]][[continuous_cov_name]])
      norm_mat <- sample_norm_mat_list[[s]][top4_genes, , drop = FALSE]
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
    
    df_spots <- do.call(rbind, spot_data_list)
    
    # Partition gradient into 10 equal-width bins
    g_range <- range(df_spots$gradient, na.rm = TRUE)
    breaks  <- seq(g_range[1], g_range[2], length.out = 11)
    df_spots$bin_idx <- cut(df_spots$gradient, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    midpoints <- (breaks[-1] + breaks[-length(breaks)]) / 2
    df_spots$bin_x <- midpoints[df_spots$bin_idx]
    
    # Calculate Mean and Standard Error (SE) per gene, active condition, and bin
    bin_summary <- df_spots %>%
      filter(!is.na(bin_x)) %>%
      group_by(gene, condition, bin_x) %>%
      summarise(
        mean_expr = mean(expression, na.rm = TRUE),
        se_expr   = sd(expression, na.rm = TRUE) / sqrt(n()),
        n_spots   = n(),
        .groups   = "drop"
      ) %>%
      mutate(
        ymin = mean_expr - se_expr,
        ymax = mean_expr + se_expr,
        gene = factor(gene, levels = top4_genes),
        condition = factor(condition, levels = active_conds)
      )
    
    # Dodge error bars only if multiple conditions are present
    bin_width <- breaks[2] - breaks[1]
    pd <- if (length(active_conds) > 1) position_dodge(width = bin_width * 0.25) else position_identity()
    
    p_grad <- ggplot(bin_summary, aes(x = bin_x, y = mean_expr, color = condition, group = condition)) +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = bin_width * 0.3, position = pd, alpha = 0.75, linewidth = 0.6) +
      geom_line(position = pd, linewidth = 1.0) +
      geom_point(position = pd, size = 2.4) +
      scale_color_manual(values = cond_palette) +
      facet_wrap(~gene, scales = "free_y", ncol = 2) +
      theme_boxed +
      labs(
        title = cd,
        x = "Motor Neuron Density",
        y = "Normalized Expression (log2)",
        color = "Condition"
      )
    
    grad_file <- file.path(gradient_dir, paste0("gradient_profile_", clean_cd, ".png"))
    ggsave(grad_file, plot = p_grad, width = 9.5, height = 7.5, dpi = 300)
    cat("    Gradient profile plot saved to:", grad_file, "\n")

  } else {
    # -------------------------------------------------------------------------
    # MODE B: Categorical Niche -> Sample x Active Compartment Z-Score Heatmap
    # -------------------------------------------------------------------------
    top_sig <- df_c %>% arrange(padj_tessera, desc(abs(estimate))) %>% head(30)
    if (nrow(top_sig) == 0) next
    
    top_sig     <- top_sig %>% arrange(estimate)
    gene_levels <- top_sig$gene
    
    all_levels    <- unique(as.character(sample_covs[[target_cov]]))
    active_levels <- all_levels[sapply(all_levels, function(lev) {
      any(grepl(paste0(target_cov, lev), active_coefs, ignore.case = TRUE))
    })]
    if (length(active_levels) == 0) active_levels <- all_levels
    
    cat(sprintf("    Generating categorical heatmap for conditions [%s] across levels [%s]...\n", 
                paste(active_conds, collapse = ", "), paste(active_levels, collapse = ", ")))
    
    cell_list <- list()
    for (s in active_samples$sample) {
      s_cond   <- active_samples$condition[active_samples$sample == s]
      cov_s    <- as.character(TESSERA_data$covariates_list[[s]][[target_cov]])
      norm_mat <- sample_norm_mat_list[[s]][gene_levels, , drop = FALSE]
      
      for (lev in active_levels) {
        spots_idx <- which(cov_s == lev)
        if (length(spots_idx) == 0) next
        
        mean_lev <- Matrix::rowMeans(norm_mat[, spots_idx, drop = FALSE])
        label_suffix <- if (length(active_levels) > 1) {
          paste(" |", ifelse(tolower(lev) == "true", "Soma", ifelse(tolower(lev) == "false", "Halo", lev)))
        } else ""
        
        cell_list[[paste0(s, "_", lev)]] <- data.frame(
          gene       = gene_levels,
          sample     = s,
          condition  = s_cond,
          row_label  = paste0("[", s_cond, "] ", s, label_suffix),
          expression = log2(mean_lev + 1),
          stringsAsFactors = FALSE
        )
      }
    }
    
    df_raw <- do.call(rbind, cell_list)
    df_heat <- df_raw %>%
      group_by(gene) %>%
      mutate(score = if (sd(expression) == 0 || is.na(sd(expression))) 0 else (expression - mean(expression)) / sd(expression)) %>%
      ungroup() %>%
      mutate(gene = factor(gene, levels = gene_levels))
    
    row_order <- rev(unique(df_heat$row_label))
    df_heat$row_label <- factor(df_heat$row_label, levels = row_order)
    
    p_heat <- ggplot(df_heat, aes(x = gene, y = row_label, fill = score)) +
          geom_tile(color = "white", linewidth = 0.6) +
          scale_fill_gradient2(
            low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0, 
            name = "Z-Score",
            # Widens the colorbar and centers the title to give tick numbers room
            guide = guide_colorbar(
              barwidth = unit(5, "cm"),
              barheight = unit(0.4, "cm"),
              title.position = "top",
              title.hjust = 0.5
            )
          ) +
          theme_boxed +
          theme(
            axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold.italic", size = 10),
            axis.text.y  = element_text(face = "bold", size = 10),
            legend.text  = element_text(size = 8),                          # Smaller tick numbers
            legend.title = element_text(face = "bold", size = 9, hjust = 0.5) # Centered legend title
          ) +
          labs(
            title = cd, 
            x = NULL,        # Removed x-axis title
            y = NULL     # Simplified y-axis title
          )
    
    plot_height <- max(2.5, 0.45 * length(unique(df_heat$row_label)) + 1.8)
    heat_file   <- file.path(heat_dir, paste0("heatmap_", clean_cd, ".png"))
    ggsave(heat_file, plot = p_heat, width = 10, height = plot_height, dpi = 300)
    cat("    Heatmap saved to:", heat_file, "\n")
  }
}

cat(sprintf("[%s] TESSERA inference and plotting completed successfully!\n", Sys.time()))


