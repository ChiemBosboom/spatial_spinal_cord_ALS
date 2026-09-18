if (!is.null(snakemake@log[[1]])) {
  log_file <- file(snakemake@log[[1]], open = "wt")
  sink(log_file)
  sink(log_file, type = "message")
}

library(Seurat)
library(dplyr)
library(Matrix)

# 1. Load Data and Parse Configuration
sn <- readRDS(snakemake@input[["seurat_rds"]])
meta_df <- sn@meta.data

celltype_column        <- snakemake@params[["celltype_column"]]
group_column           <- snakemake@params[["group_column"]]
sample_column          <- snakemake@params[["sample_column"]]
max_cells_per_celltype <- as.numeric(snakemake@params[["max_cells_per_celltype"]])

is_group_valid <- !is.null(group_column) && 
                  !is.na(group_column) && 
                  group_column != "" && 
                  group_column != "None" && 
                  group_column != "NULL"

if (is_group_valid) {
  split_cols <- c(celltype_column, group_column)
  num_groups <- length(unique(meta_df[[group_column]]))
  max_cells_per_group <- floor(max_cells_per_celltype / num_groups)
} else {
  split_cols <- c(celltype_column)
  max_cells_per_group <- max_cells_per_celltype
}

temp_id_col <- "internal_barcode_temp_id"
while (temp_id_col %in% colnames(meta_df)) {
  temp_id_col <- paste0(temp_id_col, "_unique")
}
meta_df[[temp_id_col]] <- rownames(meta_df)

# 2. Proportional Stratified Subsampling
set.seed(123)
subsampled_meta <- meta_df %>%
  group_by(across(all_of(split_cols))) %>%
  group_split() %>%
  lapply(function(df) {
    total_cells <- nrow(df)
    if (total_cells == 0 || total_cells <= max_cells_per_group) {
      return(df)
    }

    sample_split  <- split(df, df[[sample_column]], drop = TRUE)
    sample_counts <- sapply(sample_split, nrow)

    target_n <- round((sample_counts / total_cells) * max_cells_per_group)
    target_n <- pmax(1, target_n)
    target_n <- pmin(target_n, sample_counts)

    sampled_list <- lapply(seq_along(sample_split), function(i) {
      s_df <- sample_split[[i]]
      n_to_sample <- target_n[i]
      if (n_to_sample >= nrow(s_df)) {
        return(s_df)
      }
      sampled_rows <- sample(seq_len(nrow(s_df)), size = n_to_sample, replace = FALSE)
      return(s_df[sampled_rows, , drop = FALSE])
    })

    sampled_df <- bind_rows(sampled_list)
    if (nrow(sampled_df) > max_cells_per_group) {
      final_indices <- sample(seq_len(nrow(sampled_df)), size = max_cells_per_group, replace = FALSE)
      sampled_df <- sampled_df[final_indices, , drop = FALSE]
    }
    return(sampled_df)
  }) %>%
  bind_rows()

sn_subsampled <- sn[, subsampled_meta[[temp_id_col]]]

# 3. Cell Count Summary Generation
if (is_group_valid) {
  summary_matrix <- table(
    sn_subsampled@meta.data[[celltype_column]], 
    sn_subsampled@meta.data[[group_column]]
  )
  summary_df <- as.data.frame.matrix(summary_matrix)
  summary_df$total_cells <- rowSums(summary_df)
  summary_df <- data.frame(celltype = rownames(summary_df), summary_df, check.names = FALSE)
} else {
  summary_vector <- table(sn_subsampled@meta.data[[celltype_column]])
  summary_df <- data.frame(
    celltype = names(summary_vector), 
    total_cells = as.numeric(summary_vector), 
    check.names = FALSE
  )
}

write.csv(summary_df, file = snakemake@output[["summary"]], row.names = FALSE)

# 4. Export Intermediate Data
if (inherits(sn_subsampled[["RNA"]], "Assay5")) {
  sn_subsampled <- JoinLayers(sn_subsampled, assay = "RNA")
}

counts_matrix <- sn_subsampled[["RNA"]]$counts
writeMM(counts_matrix, file = snakemake@output[["mtx"]])

writeLines(rownames(sn_subsampled), con = snakemake@output[["features"]])
writeLines(colnames(sn_subsampled), con = snakemake@output[["barcodes"]])

write.csv(sn_subsampled@meta.data, file = snakemake@output[["metadata"]], row.names = TRUE)