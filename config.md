# Pipeline Configuration Guide (`config.yaml`)

This document provides a detailed reference for all configuration parameters used in the Visium HD Motor Neuron Spatial Pipeline. The configuration file governs the entire execution flow—from single-nucleus reference sub-sampling and GPU-accelerated deconvolution to geometric neuron segmentation and spatial differential expression via TESSERA.

---

## Quick Navigation
1. [Global Storage & Input Datasets](#1-global-storage--input-datasets)
2. [Stage 01: Single-Cell Reference Prep & Model Training](#2-stage-01-single-cell-reference-prep--model-training)
3. [Stage 02: Visium HD Spatial Deconvolution (cell2location)](#3-stage-02-visium-hd-spatial-deconvolution-cell2location)
4. [Stage 03: Motor Neuron Segmentation & Halo Gradients](#4-stage-03-motor-neuron-segmentation--halo-gradients)
5. [Stage 04: Cell Type Abundance Comparisons](#5-stage-04-cell-type-abundance-comparisons)
6. [Stage 05: Spatial Differential Expression (TESSERA)](#6-stage-05-spatial-differential-expression-tessera)
7. [Parameter Summary Matrix](#7-parameter-summary-matrix)

---

## 1. Global Storage & Input Datasets

```yaml
output_dir: "FINAL"
seurat_rds: "/path/to/SpinalCord_SingleNucleus_v2.rds"
gm_table_path: "/path/to/Barcodes_GreyMatter.csv"

visium_samples:
  UMC-CD-x030-s:
    path: "/path/to/spaceranger/outs/binned_outputs/square_016um/"
    condition: "ALS"
```

### `output_dir`
* **Type:** `string` (Path)
* **Description:** Base directory where all rule outputs, intermediate files, trained PyTorch/Pyro models, tables, and visualization galleries will be written.
* **Best Practice:** Use a project-specific identifier or semantic version tag (e.g., `results_v1_16um`). Avoid running multiple concurrent Snakemake instances pointing to the same `output_dir`.

### `seurat_rds`
* **Type:** `string` (Filepath)
* **Description:** Path to an annotated Seurat object (`.rds`) containing single-nucleus or single-cell RNA-seq reference data.
* **Requirements:**
  * Must contain raw counts (un-normalized integer UMI matrix) in `obj[["RNA"]]$counts` (Seurat v5) or `obj@assays$RNA@counts` (Seurat v3/v4).
  * Cell type annotations and batch metadata must be present in `meta.data`.

### `gm_table_path`
* **Type:** `string` (Filepath)
* **Description:** CSV table containing spot-level anatomical classifications (identifying spots within the spinal cord grey matter).
* **Format:** Must contain at least two columns:
  * `sample_id`: Matches keys in `visium_samples`.
  * `Barcode`: Visium HD spot coordinate barcode (e.g., `s_016um_00120_00045-1`).
* **Biological Context:** Ventral horn motor neurons reside exclusively in the grey matter. Applying this mask eliminates edge artifacts in the white matter tracts and prevents spurious segmentation calls.

### `visium_samples`
* **Type:** `mapping` (Key-Value Dictionary)
* **Description:** Dictionary defining all Visium HD capture areas to process.
  * **Key:** Unique sample identifier (e.g., `UMC-CD-x030-s`).
  * `path`: Absolute path to the Space Ranger binned output directory containing:
    * `filtered_feature_bc_matrix.h5`
    * `spatial/tissue_positions.parquet`
    * `spatial/scalefactors_json.json`
    * `spatial/tissue_hires_image.png`
  * `condition`: Experimental cohort label (e.g., `ALS`, `CTRL`). Used by statistical rules to construct design matrices and linear contrasts.

---

## 2. Stage 01: Single-Cell Reference Prep & Model Training

This stage balances single-nucleus reference cell types and trains a Negative Binomial regression model (`cell2location.models.RegressionModel`) to infer cluster-specific expression signatures ($\mu_{g,c}$).

```yaml
celltype_column: "motor_neuron_predicted"
group_column: "disease_family"
sample_column: "orig.ident"
max_cells_per_celltype: 2000

cell_count_cutoff: 5
cell_percentage_cutoff2: 0.05
nonz_mean_cutoff: 1.2

max_epochs: 700
batch_size: 16000
```

### Reference Subsampling

| Parameter | Type | Default | Rationale & Guidance |
| :--- | :--- | :--- | :--- |
| `celltype_column` | `string` | `"motor_neuron_predicted"` | Column in `seurat_rds@meta.data` containing cell type classifications. Must accurately identify the target motor neuron population. |
| `group_column` | `string` \| `null` | `"disease_family"` | Metadata column used for balanced stratified subsampling across biological batches or donor cohorts. Set to `null` to disable cohort balancing. |
| `sample_column` | `string` | `"orig.ident"` | Metadata column denoting biological library/batch. Cell2location uses this to model technical batch effects ($\epsilon_{g,e}$). |
| `max_cells_per_celltype` | `integer` | `2000` | Upper limit of cells to sample per cell type. Ensures abundant populations (e.g., oligodendrocytes) do not dominate GPU memory or bias regression fitting over rare populations (e.g., motor neurons). |

### Reference Gene Filtering

Cell2location infers average gene expression per cluster. Including unexpressed or ultra-sparse genes introduces stochastic noise into deconvolution.

* **`cell_count_cutoff` (`integer`, default: `5`):** Minimum raw count a gene must have across all cells to be considered detected.
* **`cell_percentage_cutoff2` (`float`, default: `0.05`):** Gene must be expressed in at least this proportion of cells (5%) in **at least one** cell type cluster. Retains cluster-specific markers while purging non-specific background noise.
* **`nonz_mean_cutoff` (`float`, default: `1.2`):** Minimum average count among non-zero cells. Prevents genes with solitary, low-confidence counts from inflating signature matrices.

### Model Training Parameters

* **`max_epochs` (`integer`, default: `700`):** Total training iterations. Convergence should be confirmed via the generated `training_history.png` (ELBO plateau).
* **`batch_size` (`integer`, default: `16000`):** Number of cells per mini-batch. Scaled for high-memory GPUs (e.g., NVIDIA A100 80GB). Decrease to `4096` or `2048` if training on 16GB–24GB GPUs (V100/RTX 3090).

---

## 3. Stage 02: Visium HD Spatial Deconvolution (cell2location)

This stage fits cell type abundance per bin ($w_{s,c}$) on 16 µm Visium HD grids using `cell2location.models.Cell2location`.

```yaml
umi_min: 50
gene_min: 25

max_epochs_spatial: 500
batch_size_spatial: 16000
N_cells_per_location: 1
detection_alpha: 20
```

### Spot QC Thresholds

| Parameter | Type | Default | Rationale & Guidance |
| :--- | :--- | :--- | :--- |
| `umi_min` | `integer` | `50` | Minimum UMI counts for a 16 µm bin to be included in deconvolution. Visium HD bins have significantly lower counts than legacy 55 µm Visium spots. |
| `gene_min` | `integer` | `25` | Minimum unique genes detected per bin. Eliminates bins falling on acellular tissue tears or empty glass. |

> [!TIP]
> If transitioning from **16 µm** to **8 µm** bins, decrease `umi_min` to ~`15–20` and `gene_min` to ~`10–15` to avoid dropping true cellular bins. For **2 µm** bins, specialized subcellular pipelines should be considered.

### Spatial Hyperparameters

* **`N_cells_per_location` (`integer`, default: `1`):**
  * **Mathematical Meaning:** Sets the prior mean for the total number of cells expected per spatial location ($m_s \sim \text{Gamma}(N, \dots)$).
  * **Visium HD Context:** In standard 55 µm Visium, this is typically set to `10–30`. Because a 16 µm bin approximates the cross-sectional footprint of a single eukaryotic cell body, `N_cells_per_location: 1` is mathematically and biologically appropriate.
* **`detection_alpha` (`float`, default: `20`):**
  * **Mathematical Meaning:** Regularization hyperparameter governing the variance of spot-specific sensitivity effects ($y_s$).
  * **Guidance:** A value of `20` assumes moderate technical efficiency variation across bins without allowing sensitivity differences to absorb true biological abundance shifts.
* **`max_epochs_spatial` (`integer`, default: `500`):** Training epochs for spatial mapping. Check `spatial_training_history.png` to ensure the ELBO has stabilized.
* **`batch_size_spatial` (`integer`, default: `16000`):** Spatial spots processed per optimization step. Reduce if encountering CUDA Out-of-Memory (OOM) errors.

---

## 4. Stage 03: Motor Neuron Segmentation & Halo Gradients

This stage segments individual motor neuron somas from continuous deconvolution probabilities and builds a spatial microenvironment density field.

```yaml
target_col: "Motor Neurons"
bin_size_um: 16.0
core_prob_threshold: 0.5
min_core_bins: 5
relaxed_prob_threshold: 0.25
max_expansion_um: 25.0
marker_genes:
  - "CHAT"
  - "SLC5A7"
min_marker_counts: 1
require_grey_matter: true
decay_um: 100.0
```

### Segmentation Algorithm Architecture

```
[Cell2location Abundance]
         │
         ▼
[Core Seeding] ────────► Score >= core_prob_threshold (0.5)
         │               DBSCAN clustering (eps = 1.5 bins, min_bins = 5)
         ▼
[Soma Expansion] ──────► Score >= relaxed_prob_threshold (0.25)
         │               Dijkstra geodesic expansion (dist <= max_expansion_um)
         ▼
[Dual Verification] ───► Total CHAT/SLC5A7 counts >= min_marker_counts (1)
         │               Inside Grey Matter mask (require_grey_matter = true)
         ▼
[Halo Density Field] ──► Exponential distance decay to nearest neuron:
                         exp(-distance_um / decay_um)
```

### Parameter Details

#### `target_col`
* **Type:** `string`
* **Description:** Cell type factor name produced by cell2location (derived from the reference AnnData). Must match a column in `q05_cell_abundance_w_sf`.

#### `core_prob_threshold` & `min_core_bins`
* **Defaults:** `0.5`, `5`
* **Rationale:** A candidate motor neuron seed must contain at least `5` contiguous bins each exhibiting a 5th-percentile abundance score $\ge 0.5$. Spinal motor neurons are among the largest cells in the mammalian nervous system ($40–80\ \mu\text{m}$ diameter), so a true soma inevitably covers multiple 16 µm bins.

#### `relaxed_prob_threshold` & `max_expansion_um`
* **Defaults:** `0.25`, `25.0`
* **Rationale:** Motor neuron boundaries taper outward, reducing peripheral deconvolution scores. Starting from confirmed core bins, Dijkstra’s shortest path algorithm traverses adjacent bins with abundance $\ge 0.25$ up to a Euclidean geodesic distance of $25\ \mu\text{m}$ from the core centroid.

#### `marker_genes` & `min_marker_counts`
* **Defaults:** `["CHAT", "SLC5A7"]`, `1`
* **Rationale:** Deconvolution models can occasionally produce false-positive motor neuron calls in high-density interneuron fields. To ensure fidelity, segmented candidate somas must express at least `1` raw UMI of definitive cholinergic machinery (*Choline Acetyltransferase* or the high-affinity choline transporter *SLC5A7*).

#### `require_grey_matter`
* **Type:** `boolean` (`true` / `false`)
* **Description:** Enforces that segmented soma bins reside within annotated grey matter boundaries. Removes edge-effect deconvolution false positives in the lateral or ventral white matter columns.

#### `decay_um`
* **Type:** `float` (Microns)
* **Default:** `100.0`
* **Mathematical Definition:**
  $$\text{motor\_neuron\_density}_s = \exp\left(-\frac{\text{dist}(s, \mathcal{MN})}{\text{decay\_um}}\right)$$
  where $\text{dist}(s, \mathcal{MN})$ is the physical Euclidean distance (in $\mu\text{m}$) from bin $s$ to the boundary of the nearest verified motor neuron soma.
* **Biological Context:** At $100\ \mu\text{m}$, density equals $e^{-1} \approx 0.368$; at $200\ \mu\text{m}$, it drops to $e^{-2} \approx 0.135$. This models paracrine signaling, neuroinflammatory gradients, and glial activation surrounding degenerating motor neurons.

---

## 5. Stage 04: Cell Type Abundance Comparisons

This stage computes sample-level and condition-level compositional proportions across defined tissue compartments using `scripts/compare_cells.R`.

```yaml
bin_comparisons:
  is_grey_matter:
    label: "Grey Matter"
    filter: "is_grey_matter"
  motor_neuron_density:
    label: "MN Microenvironment"
    filter: "motor_neuron_density > 0.25"
```

### Defining Compartment Comparisons
Each entry in `bin_comparisons` generates a dedicated set of abundance strip-plots, stacked bar charts, and summary statistics:
* **Key:** Machine-readable identifier for the comparison (used in file naming).
* `label`: Clean string displayed on publication plot headers.
* `filter`: A valid R/`dplyr` logical expression evaluated against spot metadata (`adata_vis.obs`).

### Example Filters
```yaml
# Broad grey matter background
filter: "is_grey_matter"

# Immediate perineuronal microenvironment (~140 um halo)
filter: "motor_neuron_density > 0.25"

# Intimate contact soma/juxtaneuronal zone
filter: "is_neuron_expanded == TRUE"

# Distal grey matter controls (far from any motor neuron)
filter: "is_grey_matter & motor_neuron_density < 0.05"
```

---

## 6. Stage 05: Spatial Differential Expression (TESSERA)

TESSERA fits a generalized linear mixed model (GLMM) with a Leroux conditional autoregressive (CAR) random effect to account for spatial autocorrelation across spots while testing disease covariates.

```yaml
tessera_d_thresh: 23.0
tessera_min_pct_spots: 0.01
tessera_min_nonz_mean_counts: 1.1
tessera_fdr_threshold: 0.05

tessera_analyses:
  mn_gradient:
    label: "MN Density Gradient"
    filter: "motor_neuron_density > 0.01 & motor_neuron_density < 1"
    design_formula: "~ 0 + condition + condition:motor_neuron_density"
    contrasts:
      slope_ALS_vs_CTRL: "conditionALS:motor_neuron_density - conditionCTRL:motor_neuron_density"
      slope_ALS: "conditionALS:motor_neuron_density"
      slope_CTRL: "conditionCTRL:motor_neuron_density"

  mn_microenvironment:
    label: "MN Microenvironment"
    filter: "motor_neuron_density > 0.25"
    design_formula: "~ 0 + condition"
    contrasts:
      ALS_vs_CTRL: "conditionALS - conditionCTRL"
```

### Global TESSERA Parameters

* **`tessera_d_thresh` (`float`, default: `23.0`):**
  * **Mathematical Meaning:** Maximum Euclidean distance (in image full-resolution pixel space) connecting two bins as spatial neighbors in the adjacency graph $W$.
  * **Calibration:** In Visium HD 16 µm grids, bin centers are regularized. A threshold of `23.0` connects immediate orthogonal and diagonal neighbors (first-order spatial lag) without creating an overly dense adjacency graph.
* **`tessera_min_pct_spots` (`float`, default: `0.01`):** Minimum detection fraction (1% of analyzed spots). Genes expressed in fewer spots are dropped to prevent GLMM non-convergence.
* **`tessera_min_nonz_mean_counts` (`float`, default: `1.1`):** Average raw count within spots expressing the gene ($>0$). Eliminates ultra-low expression genes characterized by single UMIs.
* **`tessera_fdr_threshold` (`float`, default: `0.05`):** Significance threshold applied to empirical null-adjusted $p$-values.

### Statistical Models & Linear Contrasts

The `tessera_analyses` block defines modular differential expression runs.

#### Analysis 1: `mn_gradient` (Continuous Slope Test)
* **Goal:** Test whether the rate of gene expression change as a function of distance to motor neurons differs between ALS and CTRL.
* **Filter:** `motor_neuron_density > 0.01 & motor_neuron_density < 1` (excludes the soma core itself and distal tissue).
* **Design Formula:** `~ 0 + condition + condition:motor_neuron_density`
  * `conditionALS` and `conditionCTRL`: Condition-specific intercepts.
  * `conditionALS:motor_neuron_density`: Slope of expression with respect to motor neuron density in ALS.
  * `conditionCTRL:motor_neuron_density`: Slope of expression with respect to motor neuron density in CTRL.
* **Contrast `slope_ALS_vs_CTRL`:**
  $$\beta_{\text{ALS:density}} - \beta_{\text{CTRL:density}} = 0$$
  Identifies genes whose spatial proximity gradients are significantly steepened or flattened in ALS.

#### Analysis 2: `mn_microenvironment` (Discrete Niche Test)
* **Goal:** Test differential expression exclusively within the perineuronal niche.
* **Filter:** `motor_neuron_density > 0.25` (bins within $\sim 140\ \mu\text{m}$ of a motor neuron).
* **Design Formula:** `~ 0 + condition`
* **Contrast `ALS_vs_CTRL`:**
  $$\beta_{\text{ALS}} - \beta_{\text{CTRL}} = 0$$
  Calculates standard log2 fold change within the immediate microenvironment while controlling for spatial autocorrelation.

---

## 7. Parameter Summary Matrix

| Parameter | Recommended Default | Safe Range | Primary Impact |
| :--- | :--- | :--- | :--- |
| `max_cells_per_celltype` | `2000` | `500 – 5000` | Reference balance and GPU training duration |
| `cell_count_cutoff` | `5` | `1 – 20` | Reference gene sparsity filter |
| `nonz_mean_cutoff` | `1.2` | `1.05 – 1.5` | Removes low-count gene noise from signatures |
| `umi_min` (16 µm) | `50` | `25 – 100` | Visium HD spot quality control |
| `gene_min` (16 µm) | `25` | `15 – 50` | Visium HD spot quality control |
| `N_cells_per_location` | `1` | `1 – 2` | Cell2location prior for 16 µm bins |
| `detection_alpha` | `20` | `10 – 100` | Spot sensitivity regularization |
| `core_prob_threshold` | `0.5` | `0.3 – 0.7` | Stringency of motor neuron soma seeds |
| `min_core_bins` | `5` | `3 – 10` | Minimum physical footprint of a neuron core |
| `max_expansion_um` | `25.0` | `10.0 – 40.0` | Maximum radius of soma expansion |
| `decay_um` | `100.0` | `50.0 – 250.0` | Spatial extent of the perineuronal density field |
| `tessera_d_thresh` | `23.0` | `18.0 – 30.0` | Spatial graph neighborhood connectivity |
| `tessera_min_pct_spots`| `0.01` | `0.005 – 0.05` | Filter for stable GLMM convergence |
| `tessera_fdr_threshold`| `0.05` | `0.01 – 0.10` | False discovery rate cutoff |
