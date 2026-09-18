#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Setup Logging & Environment
# -----------------------------------------------------------------------------
log_file <- snakemake@log[[1]]
if (!is.null(log_file) && log_file != "") {
  log_con <- file(log_file, open = "wt")
  sink(log_con)
  sink(log_con, type = "message")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(arrow)
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(S4Vectors)
  library(Matrix)
  library(MatrixGenerics)
  library(TESSERA)
  library(future)
  library(future.apply)
  library(parallel)
})

# Prevent BLAS thread oversubscription / deadlocks when forking worker processes
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

# Inputs & Outputs
metadata_path    <- snakemake@input[["metadata"]]
visium_samples   <- snakemake@params[["visium_samples"]]
out_tessera_data <- snakemake@output[["tessera_data"]]
out_tessera_fits <- snakemake@output[["tessera_fits"]]
out_perf_summary <- snakemake@output[["perf_summary"]]

# Parameters
umi_min             <- as.numeric(snakemake@params[["umi_min"]])
gene_min            <- as.numeric(snakemake@params[["gene_min"]])
d_thresh            <- as.numeric(snakemake@params[["d_thresh"]])
min_pct_spots       <- as.numeric(snakemake@params[["min_pct_spots"]])
min_nonzero_mean    <- as.numeric(snakemake@params[["min_nonz_mean_counts"]])

filter_expr         <- snakemake@params[["filter_expr"]]
design_formula      <- snakemake@params[["design_formula"]]

# TESSERA Default Modeling Configurations
model_type   <- "Leroux"
em_iters     <- 200
opt_iters    <- 5
em_min_iters <- 30
em_tol       <- 1e-3
em_stopping  <- "rel_loglike"
beta_init    <- "glm"
gamma_init   <- "moran"
tau2_init    <- "var"

cat(sprintf("[%s] Starting TESSERA model fitting pipeline...\n", Sys.time()))
cat(sprintf("D_THRESH: %.1f | Design Formula: %s\n", d_thresh, design_formula))

# -----------------------------------------------------------------------------
# 1. Load and QC a Single Visium HD Sample
# -----------------------------------------------------------------------------
load_visium_hd_sample <- function(sample_id, sample_info, umi_min, gene_min) {
  base_path <- sample_info$path
  condition <- sample_info$condition

  h5_path      <- file.path(base_path, "filtered_feature_bc_matrix.h5")
  parquet_path <- file.path(base_path, "spatial/tissue_positions.parquet")

  if (!file.exists(h5_path)) stop(paste("Count matrix not found at:", h5_path))
  if (!file.exists(parquet_path)) stop(paste("Tissue positions not found at:", parquet_path))

  cat(sprintf("  Loading sample '%s' (Condition: %s)...\n", sample_id, condition))

  counts <- Seurat::Read10X_h5(h5_path)
  coords <- as.data.frame(arrow::read_parquet(parquet_path))
  rownames(coords) <- coords$barcode

  common_bc <- intersect(colnames(counts), coords$barcode)
  counts    <- counts[, common_bc, drop = FALSE]
  coords    <- coords[common_bc, , drop = FALSE]

  n_umi   <- Matrix::colSums(counts)
  n_gene  <- Matrix::colSums(counts > 0)
  qc_keep <- (n_umi >= umi_min) & (n_gene >= gene_min)

  counts <- counts[, qc_keep, drop = FALSE]
  coords <- coords[qc_keep, , drop = FALSE]

  cat(sprintf("    Retained %d / %d spots passing QC\n", ncol(counts), length(common_bc)))

  unique_bc <- paste0(sample_id, "_", colnames(counts))
  colnames(counts) <- unique_bc
  rownames(coords) <- unique_bc

  spatial_coords <- as.matrix(coords[, c("pxl_col_in_fullres", "pxl_row_in_fullres")])
  colnames(spatial_coords) <- c("pxl_col_in_fullres", "pxl_row_in_fullres")
  rownames(spatial_coords) <- unique_bc

  col_data <- DataFrame(
    barcode   = colnames(counts),
    sample_id = sample_id,
    condition = condition,
    nUMI      = Matrix::colSums(counts),
    nGene     = Matrix::colSums(counts > 0),
    row.names = unique_bc
  )

  spe <- SpatialExperiment(
    assays        = list(counts = counts),
    colData       = col_data,
    spatialCoords = spatial_coords,
    sample_id     = sample_id
  )

  return(spe)
}

# -----------------------------------------------------------------------------
# 2. Combine Samples into Joint SpatialExperiment
# -----------------------------------------------------------------------------
spe_list <- lapply(names(visium_samples), function(s_name) {
  load_visium_hd_sample(s_name, visium_samples[[s_name]], umi_min, gene_min)
})
names(spe_list) <- names(visium_samples)

cat("Combining samples into joint SpatialExperiment...\n")
sample_gene_lists <- lapply(spe_list, rownames)
all_union_genes   <- length(Reduce(union, sample_gene_lists))
common_genes      <- Reduce(intersect, sample_gene_lists)

cat(sprintf("  Sample gene intersection: retained %d common genes (filtered out %d genes not shared by all samples).\n",
            length(common_genes), all_union_genes - length(common_genes)))

spe_list  <- lapply(spe_list, function(spe) spe[common_genes, ])
joint_spe <- do.call(cbind, spe_list)
cat(sprintf("Joint dimensions: %d genes x %d spots\n", nrow(joint_spe), ncol(joint_spe)))

# -----------------------------------------------------------------------------
# 3. Integrate Spot-Level Metadata
# -----------------------------------------------------------------------------
cat("Reading spot metadata from:", metadata_path, "\n")
spot_metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
rownames(spot_metadata) <- spot_metadata$spot_id

common_spots <- intersect(colnames(joint_spe), rownames(spot_metadata))
if (length(common_spots) == 0) {
  stop("FATAL: No overlapping barcodes between Visium counts and spot metadata!")
}
cat(sprintf("Matched %d spots between count matrix and metadata.\n", length(common_spots)))

joint_spe    <- joint_spe[, common_spots]
aligned_meta <- spot_metadata[common_spots, , drop = FALSE]

for (col in colnames(aligned_meta)) {
  colData(joint_spe)[[col]] <- aligned_meta[[col]]
}

# -----------------------------------------------------------------------------
# 4. Subset Spots via Analysis Filter Expression
# -----------------------------------------------------------------------------
if (!is.null(filter_expr) && filter_expr != "" && filter_expr != "all") {
  cat(sprintf("Filtering spots using: '%s'...\n", filter_expr))
  df_meta <- as.data.frame(colData(joint_spe))
  keep_spots_idx <- eval(parse(text = filter_expr), envir = df_meta)
  keep_spots_idx[is.na(keep_spots_idx)] <- FALSE

  joint_spe <- joint_spe[, keep_spots_idx]
  cat(sprintf("Spots remaining after filter: %d\n", ncol(joint_spe)))

  if (ncol(joint_spe) == 0) {
    stop("FATAL: Zero spots remained after applying the metadata filter.")
  }
}

# -----------------------------------------------------------------------------
# 5. Covariate Formatting & Missing Value Validation
# -----------------------------------------------------------------------------
cat("Preparing metadata variables for design matrix...\n")

joint_spe$nested_id <- "1"
for (cond in unique(joint_spe$condition)) {
  samps_in_cond <- unique(joint_spe$sample_id[joint_spe$condition == cond])
  for (idx in seq_along(samps_in_cond)) {
    s_id <- samps_in_cond[idx]
    joint_spe$nested_id[joint_spe$sample_id == s_id] <- as.character(idx)
  }
}
joint_spe$nested_id <- factor(joint_spe$nested_id)

formula_obj  <- stats::as.formula(design_formula)
formula_vars <- all.vars(formula_obj)

df_vars <- as.data.frame(colData(joint_spe)[, formula_vars, drop = FALSE])
complete_spots <- stats::complete.cases(df_vars)
if (any(!complete_spots)) {
  cat(sprintf("  Removing %d spots with NA values in formula variables.\n", sum(!complete_spots)))
  joint_spe <- joint_spe[, complete_spots]
  df_vars   <- df_vars[complete_spots, , drop = FALSE]
}

for (v in formula_vars) {
  val <- df_vars[[v]]
  if (is.logical(val) || (is.character(val) && all(unique(na.omit(val)) %in% c("True", "False", "TRUE", "FALSE")))) {
    colData(joint_spe)[[v]] <- factor(as.character(val))
  } else if (is.character(val) || is.factor(val)) {
    colData(joint_spe)[[v]] <- droplevels(factor(val))
  } else if (is.numeric(val)) {
    colData(joint_spe)[[v]] <- as.numeric(val)
  }
}

# -----------------------------------------------------------------------------
# 6. Gene Filtering
# -----------------------------------------------------------------------------
cat("\n--- Filtering Genes ---\n")
initial_genes <- nrow(joint_spe)
cat(sprintf("Initial genes: %d\n", initial_genes))

# Step 6a: Mitochondrial Gene Removal
mt_mask   <- grepl("^(?i)mt-", rownames(joint_spe))
n_mt      <- sum(mt_mask)
joint_spe <- joint_spe[!mt_mask, ]
cat(sprintf("  [1/2] Mitochondrial filter: removed %d genes (Remaining: %d)\n", 
            n_mt, nrow(joint_spe)))

counts_mat <- SingleCellExperiment::counts(joint_spe)

# Step 6b: Low-Expression and Low-Prevalence Filtering (Non-zero Mean & Spot Fraction)
n_nonzero_spots     <- Matrix::rowSums(counts_mat > 0)
pct_expressed       <- n_nonzero_spots / ncol(counts_mat)
total_counts        <- Matrix::rowSums(counts_mat)

# Calculate average count only in spots where the gene is expressed (>0)
nonzero_mean_counts <- ifelse(n_nonzero_spots > 0, total_counts / n_nonzero_spots, 0)

fail_pct_spots      <- pct_expressed < min_pct_spots
fail_nonzero_mean   <- nonzero_mean_counts < min_nonzero_mean
gene_keep           <- (!fail_pct_spots) & (!fail_nonzero_mean)

cat(sprintf("  [2/2] Expression QC filters:\n"))
cat(sprintf("        - Failed spot detection rate (< %s): %d genes\n", min_pct_spots, sum(fail_pct_spots)))
cat(sprintf("        - Failed non-zero mean counts (< %g): %d genes\n", min_nonzero_mean, sum(fail_nonzero_mean)))
cat(sprintf("        - Total lowly expressed genes removed: %d (Remaining: %d)\n", 
            sum(!gene_keep), sum(gene_keep)))

joint_spe  <- joint_spe[gene_keep, ]
counts_mat <- SingleCellExperiment::counts(joint_spe)

# Step 6c: Remove Spots Without Any Counts Left
spot_totals <- Matrix::colSums(counts_mat)
valid_spots <- spot_totals > 0
if (any(!valid_spots)) {
  cat(sprintf("  Removing %d spots with 0 counts after gene filtering.\n", sum(!valid_spots)))
  joint_spe  <- joint_spe[, valid_spots]
  counts_mat <- SingleCellExperiment::counts(joint_spe)
}

# All remaining genes will be fitted
genes_to_fit <- rownames(joint_spe)
cat(sprintf("\nGene filtering complete. Total genes to fit: %d (across %d spots).\n", 
            length(genes_to_fit), ncol(joint_spe)))

# -----------------------------------------------------------------------------
# 7. Design Matrix Construction & Rank Check
# -----------------------------------------------------------------------------
df_meta <- as.data.frame(colData(joint_spe))
X_mat <- stats::model.matrix(stats::as.formula(design_formula), data = df_meta)
X_mat <- X_mat[, colSums(abs(X_mat)) > 0, drop = FALSE]

rk <- Matrix::rankMatrix(X_mat)[1]
if (rk != ncol(X_mat)) {
  stop(sprintf("Design matrix is not full rank (rank = %d, cols = %d). Check for collinearity.", rk, ncol(X_mat)))
}

# -----------------------------------------------------------------------------
# 8. Data Container Construction (TESSERA::prep_data)
# -----------------------------------------------------------------------------
cat("Executing TESSERA::prep_data...\n")
t0 <- Sys.time()

TESSERA_data <- TESSERA::prep_data(
  x          = joint_spe,
  sample_col = "sample_id",
  design_mat = X_mat,
  model_type = model_type,
  D_THRESH   = d_thresh
)

# Subset counts to all remaining QC-filtered genes
for (s_id in names(TESSERA_data)) {
  TESSERA_data[[s_id]]$counts <- TESSERA_data[[s_id]]$counts[genes_to_fit, , drop = FALSE]
}

cat(sprintf("TESSERA::prep_data finished in %.2f minutes.\n", 
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

dir.create(dirname(out_tessera_data), recursive = TRUE, showWarnings = FALSE)
saveRDS(TESSERA_data, file = out_tessera_data)

# -----------------------------------------------------------------------------
# 9. Parallel Spatial GLMM Fitting (TESSERA_lattice)
# -----------------------------------------------------------------------------
n_threads <- as.integer(snakemake@threads)
cat(sprintf("\n[%s] Fitting %d genes in parallel (%d threads)...\n", 
            Sys.time(), length(genes_to_fit), n_threads))

if (.Platform$OS.type == "unix") {
  future::plan(future::multicore, workers = n_threads)
} else {
  future::plan(future::multisession, workers = n_threads)
}
options(future.globals.maxSize = 4000 * 1024^2)

t_fit_start <- Sys.time()

fits_list <- future.apply::future_lapply(
  genes_to_fit,
  function(g_name) {
    tryCatch({
      TESSERA::TESSERA_lattice(
        TESSERAData_obj = TESSERA_data,
        gene_name       = g_name,
        model_type      = model_type,
        em_iters        = em_iters,
        opt_iters       = opt_iters,
        em_min_iters    = em_min_iters,
        em_tol          = em_tol,
        em_stopping     = em_stopping,
        beta_init       = beta_init,
        gamma_init      = gamma_init,
        tau2_init       = tau2_init,
        verbose         = FALSE
      )
    }, error = function(e) {
      warning(sprintf("Gene '%s' failed: %s", g_name, conditionMessage(e)))
      return(NULL)
    })
  },
  future.seed = TRUE,
  future.chunk.size = 1
)

future::plan(future::sequential)
names(fits_list) <- genes_to_fit

cat(sprintf("Fitting finished in %.2f minutes.\n", 
            as.numeric(difftime(Sys.time(), t_fit_start, units = "mins"))))

# -----------------------------------------------------------------------------
# 10. Collate Summaries and Save Outputs
# -----------------------------------------------------------------------------
is_converged <- !vapply(fits_list, is.null, logical(1))
n_success    <- sum(is_converged)
cat(sprintf("Successfully converged: %d / %d genes (%.1f%%).\n", 
            n_success, length(genes_to_fit), (n_success / length(genes_to_fit)) * 100))

if (n_success == 0) {
  stop("FATAL: All gene models failed to fit.")
}

TESSERA_fits <- fits_list[is_converged]

perf_list <- lapply(names(TESSERA_fits), function(g) {
  df <- TESSERA_fits[[g]]$performanceSummary
  if (!is.null(df) && nrow(df) > 0) df$gene <- g
  return(df)
})
perf_df <- do.call(rbind, perf_list)

# Fallback to an empty dataframe with column header if performance summaries are empty
if (is.null(perf_df)) {
  perf_df <- data.frame(gene = character(0))
}

dir.create(dirname(out_tessera_fits), recursive = TRUE, showWarnings = FALSE)
saveRDS(TESSERA_fits, file = out_tessera_fits)

dir.create(dirname(out_perf_summary), recursive = TRUE, showWarnings = FALSE)
write.csv(perf_df, file = out_perf_summary, row.names = FALSE)

cat(sprintf("[%s] TESSERA model fitting successfully completed!\n", Sys.time()))