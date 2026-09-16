### 1. Reference Downsampling (`downsample_seurat`)

**Overview**  
Subsamples a single-cell/nucleus Seurat object to create a balanced reference dataset for Cell2location deconvolution. It prevents overrepresented cell types from dominating model training while preserving donor diversity and balancing across experimental conditions.

**Inputs & Outputs**
* **Input:** `seurat_rds` — Path to your annotated Seurat object (`.rds`) containing raw RNA counts.
* **Output:** `cell_counts_summary.csv` — Table showing the final number of retained cells per cell type (broken down by group, if enabled).

**Configuration (`config.yaml`)**
* `celltype_column`: Metadata column containing cell-type annotations (e.g., `"motor_neuron_predicted"`).
* `group_column`: Metadata column to balance across (e.g., `"disease_family"`). Set to `null` if you do not want to stratify by group.
* `sample_column`: Metadata column for donor/sample IDs (e.g., `"orig.ident"`). Ensures downsampling pulls cells proportionally across all donors.
* `max_cells_per_celltype`: Maximum cell budget per cell type (e.g., `2000`). If a `group_column` is provided, this budget is divided evenly across groups.
