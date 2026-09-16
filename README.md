# Script: downsample_seurat.R

## Overview
`downsample_seurat.R` performs balanced, stratified subsampling on a single-cell or single-nucleus Seurat object to create a reference dataset for cell-type deconvolution (e.g., Cell2location). 

The script balances representations across cell types and experimental groups while maintaining proportional donor/sample contributions. It outputs the downsampled count matrix and annotations in standard Matrix Market format alongside a cell-count summary.

---

## Workflow Logic
1. **Grouping Assessment:** Checks if a valid `group_column` (e.g., disease condition) is supplied. If present, the target cell budget (`max_cells_per_celltype`) is divided evenly among groups.
2. **Proportional Sampling:** Within each cell-type and group subset, cells are sampled proportionally across donors/samples (`sample_column`), ensuring smaller donors contribute at least one cell where possible.
3. **Data Export:** Subsets the Seurat object, joins layers if using Seurat v5 (`Assay5`), and writes the sparse count matrix, feature list, barcode list, cell metadata, and a contingency summary table.

---

## Configurables & Parameters

### Snakemake Inputs (`snakemake@input`)
* `seurat_rds`: Path to the input `.rds` file containing a Seurat object with raw RNA counts.
* `tessera_ready`: Sentinel file indicating upstream R environment setup is complete.

### Snakemake Outputs (`snakemake@output`)
* `mtx`: Matrix Market count matrix (`.mtx`) containing raw UMI counts for subsampled cells.
* `features`: Line-delimited text file containing gene/feature names.
* `barcodes`: Line-delimited text file containing cell barcodes.
* `metadata`: CSV containing cell metadata for the subsampled cells.
* `summary`: CSV summary table reporting the final number of cells retained per cell type (stratified by group, if applicable).

### Parameters (`snakemake@params`)
| Parameter | Type | Description | Example |
| :--- | :--- | :--- | :--- |
| `celltype_column` | `character` | Metadata column identifying cell-type annotations. | `"motor_neuron_predicted"` |
| `group_column` | `character` or `NULL` | Metadata column for condition/family balancing. Set to `null` to bypass group stratification. | `"disease_family"` |
| `sample_column` | `character` | Metadata column identifying individual donors/samples for proportional representation. | `"orig.ident"` |
| `max_cells_per_celltype` | `numeric` | Upper limit of cells to retain per cell type across all samples/groups. | `2000` |
