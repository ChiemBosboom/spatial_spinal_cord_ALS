# Title

description

---

## Overview

* **`scripts/downsample_seurat.R`:** Subsamples single-nucleus reference data across cell types batches, and samples to create balanced count matrices for reference modeling.
* **`scripts/train_reference.py`:** Filters genes and trains a `cell2location` model to derive cell type-specific reference signatures from single-nucleus data.
* **`scripts/train_spatial.py`:** Fits spatial `cell2location` models to infer cell type abundance across tissue slides.
* **`scripts/segment_neurons.py`:** segments motor neurons based on `cell2location` abundance estimates, validates with cholinergic markers, and models an exponential distance-decay density field.
* **`scripts/compare_cells.R`:** Quantifies cell type proportions within defined compartments and evaluates condition-level shifts.
* **`scripts/fit_tessera.R`:** 
* **`scripts/compare_genes.R`:** 

### Pilot Study & Scalability
* **Pilot Data & Technical Limitations:** The `output_pilot` directory contains results from a 4-sample pilot cohort (**3 ALS, 1 CTRL**) run with `config.yaml` using the default parameters documented below. Most ALS-related differential expression tests yielded zero significant genes. Because the pilot contained only a single CTRL sample, nested donor terms (`condition:nested_id`) could not be included in the design formula to estimate inter-individual baseline variance. Consequently, `TESSERA` could not account for biological variation across ALS donors, which distorted the background null distribution—whereas the single-sample CTRL analysis ran without issue.
* **Future Cohorts:** In larger cohorts with multiple controls, adding `condition:nested_id` will properly control for patient-to-patient variance and enable reliable differential testing. The Snakemake pipeline readily scales to future cohorts simply by updating `config.yaml`.

---

## AI Usage

Generative artificial intelligence (**Google Gemini**) was used as an assistant to write, refactor, and optimize the Python, R, and Snakemake scripts in this repository. All outputs were manually validated.

---

## Global & Input Settings

```yaml
output_dir: "output_pilot"
seurat_rds: "/path/to/SpinalCord_SingleNucleus.rds"
gm_table_path: "/path/to/Barcodes_GreyMatter.csv"

visium_samples:
  UMC-CD-x030-s:
    path: "/path/to/spaceranger/outs/binned_outputs/square_016um/"
    condition: "ALS"
  UMC-CD-x032-s:
    path: "/path/to/spaceranger/outs/binned_outputs/square_016um/"
    condition: "ALS"
  UMC-CD-x033-s:
    path: "/path/to/spaceranger/outs/binned_outputs/square_016um/"
    condition: "CTRL"
  UMC-CD-x034-s:
    path: "/path/to/spaceranger/outs/binned_outputs/square_016um/"
    condition: "ALS"
```

### Explanation
* **`output_dir` (`string`):** Base directory where all rule outputs, intermediate files, models, tables, and figures are stored (e.g., `"output_pilot"`).
* **`seurat_rds` (`string`):** Path to the single-nucleus Seurat RDS object containing raw UMI counts and cell type labels.
* **`gm_table_path` (`string`):** Path to a CSV table containing spot-level anatomical classifications. Must contain `sample_id` and `Barcode` columns mapping spots located within the spinal cord grey matter.
* **`visium_samples` (`map`):** Dictionary specifying each Visium HD sample to process:
  * **Sample Key (e.g., `UMC-CD-x030-s`):** Unique identifier used in sample-level file naming.
  * **`path` (`string`):** Directory path to the Space Ranger binned output containing `square_016um/` files (`filtered_feature_bc_matrix.h5`, `tissue_positions.parquet`, etc.).
  * **`condition` (`string`):** Biological cohort or group label (e.g., `ALS`, `CTRL`) used to construct design matrices and contrasts.

---

## 01: Single-Nucleus Reference Training

### Configuration Parameters

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `celltype_column` | `string` | `"motor_neuron_predicted"` | Column in Seurat metadata containing cell type annotations. |
| `group_column` | `string` \| `null`| `"disease_family"` | Metadata column used for balanced subsampling across groups. Set to `null` to disable. |
| `sample_column` | `string` | `"orig.ident"` | Metadata column used for balanced subsampling across samples. |
| `max_cells_per_celltype` | `integer` | `2000` | Maximum cells sampled per cell type to speed up training. |
| `cell_count_cutoff` | `integer` | `5` | Minimum number of cells a gene must be detected in. |
| `nonz_mean_cutoff` | `float` | `1.2` | Minimum average expression of a gene within non-zero expressing cells. |
| `cell_percentage_cutoff2` | `float` | `0.05` | Genes detected in at least this fraction of cells are kept regardless of nonz_mean_cutoff. |
| `max_epochs` | `integer` | `700` | Total training epochs for the reference regression model. |
| `batch_size` | `integer` | `16000` | Batch size reference regression model. |

### Outputs (`01_reference/`)

| Output Target | Description |
| :--- | :--- |
| `tables/cell_counts_summary.csv` | Cell counts per cell type and group retained after stratified subsampling. |
| `reference.h5ad` | Downsampled reference AnnData containing exported posterior expression signatures. |
| `models/reference_model/` | Saved PyTorch/Pyro weights and parameters for the trained cell2location regression model. |
| `plots/filtering_summary.png` | Diagnostic plot showing gene inclusion. |
| `plots/training_history.png` | ELBO loss convergence curve across reference training epochs. |
| `plots/training_qc_reconstruction.png` | QC scatter plot evaluating observed vs. model-reconstructed count distributions. |
| `plots/training_qc_expression.png` | QC scatter plot comparing inferred cluster expression signatures against empirical averages. |

---

## 02: Spatial Deconvolution

### Configuration Parameters

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `umi_min` | `integer` | `50` | Minimum total UMI count per 16 µm bin. |
| `gene_min` | `integer` | `25` | Minimum unique genes detected per bin. |
| `max_epochs_spatial` | `integer` | `500` | Total training epochs for spatial mapping. |
| `batch_size_spatial` | `integer` | `16000` | Batch size for spatial mapping. |
| `N_cells_per_location` | `integer` | `1` | Prior mean for the number of cells expected per bin. |
| `detection_alpha` | `integer` | `20` | Regularization strength of per-location normalization for RNA detection sensitivity |

### Outputs (`02_deconvolution/`)

| Output Target | Description |
| :--- | :--- |
| `spatial.h5ad` | Combined spatial AnnData containing cell abundance estimates across all slides. |
| `models/spatial_model/` | Saved model weights and parameters for the trained spatial cell2location model. |
| `plots/combined_qc_umi.png` | Histogram of total UMIs across all bins with the cutoff threshold marked. |
| `plots/combined_qc_genes.png` | Histogram of detected genes across all bins with the cutoff threshold marked. |
| `plots/spatial_training_history.png` | ELBO loss convergence curve during spatial model training. |
| `plots/spatial_qc_reconstruction.png` | QC plot evaluating spatial model reconstruction accuracy. |
| `tables/cell_abundance_q05.csv` | 5th-percentile cell abundance estimates per bin for all reference cell types. |
| `tables/spatial_summary_stats.csv` | Summary cell abundance statistics across spatial bins. |
| `plots/celltypes/*_spatial_abundance.png` | Multi-slide spatial abundance overlay maps generated individually for each cell type. |

---

## 03: Motor Neuron Segmentation

### Configuration Parameters

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `target_col` | `string` | `"Motor Neurons"` | Factor name from cell2location deconvolution representing the target cell type. |
| `bin_size_um` | `float` | `16.0` | Size of bins in micrometers. |
| `core_prob_threshold` | `float` | `0.5` | Minimum abundance required to consider a bin as part of a neuron core. |
| `min_core_bins` | `integer` | `5` | Minimum direct neighbor bins passing `core_prob_threshold` required to seed a candidate neuron. |
| `relaxed_prob_threshold`| `float` | `0.25` | Lower abundance required for neuron expansion from the core. |
| `max_expansion_um` | `float` | `25.0` | Maximum physical radius (in µm) to expand outward from the core. |
| `marker_genes` | `list` | `["CHAT", "SLC5A7"]` | Canonical marker genes required to validate candidate segmented neurons. |
| `min_marker_counts` | `integer` | `1` | Minimum combined raw UMIs of marker genes required in the candidate neuron. |
| `require_grey_matter` | `boolean` | `true` | When `true`, discards candidate neurons falling outside grey matter. |
| `decay_um` | `float` | `100.0` | Distance decay rate ($\lambda$) in the exponential formula: $\text{density} = \exp(-\text{dist} / \lambda)$. Density is calculated strictly over grey matter bins (non-grey matter bins remain 0). |

### Outputs (`03_segmentation/`)

| Output Target | Description |
| :--- | :--- |
| `tables/motor_neuron_metadata_per_spot.csv` | Spot-level metadata with segmentation flags (`is_neuron_core`, `is_neuron_expanded`, `neuron_id`, `motor_neuron_density`, `is_grey_matter`). |
| `plots/motor_neurons_{sample}_bins.png` | Histology image showing all verified segmented motor neuron soma bins per sample. |
| `plots/motor_neurons_{sample}_density.png` | Continuous spatial gradient map illustrating the exponential decay density field ($0 \to 1$) per sample. |

---

## 04: Cell Type Abundance Comparisons

```yaml
bin_comparisons:
  is_grey_matter:
    label: "Grey Matter"
    filter: "is_grey_matter"
  motor_neuron_density:
    label: "MN Microenvironment"
    filter: "motor_neuron_density > 0.25"
```

### Explanation 
* **Comparison Key (e.g., `is_grey_matter`, `motor_neuron_density`):** Unique internal identifier for the comparison, used in output file names (e.g., `stats_is_grey_matter.csv`).
* **`label` (`string`):** Clean descriptive title displayed on the generated stacked bar charts and strip plots.
* **`filter` (`string`):** A valid R logical expression evaluated on the spot metadata columns:
  * `"is_grey_matter"`: Evaluates all bins within the anatomically defined grey matter.
  * `"motor_neuron_density > 0.25"`: Restricts analysis to the immediate perineuronal microenvironment (~140 µm radius around motor neurons).

### Outputs (`04_cell_comparison/`)

| Output Target | Description |
| :--- | :--- |
| `plots/cell_abundance_stacked_barchart.png` | Stacked bar chart showing condition-averaged cell type proportions across all comparisons. |
| `plots/abundance_stripplots/abundance_*.png` | Sample-level strip plots with condition mean overlays generated for each individual cell type. |
| `tables/stats_{comparison}.csv` | Condition-level summary metrics for a specific comparison. |
| `tables/stats_all_comparisons.csv` | Master summary table collating cell type abundance statistics across all configured comparisons. |

---

## 05: Spatial Differential Expression

### Global Parameters

| Parameter | Type | Default | Rationale & Guidance |
| :--- | :--- | :--- | :--- |
| `tessera_d_thresh` | `float` | `23.0` | Maximum distance defining 1st-order orthogonal and diagonal spatial neighbors. |
| `tessera_min_pct_spots` | `float` | `0.01` | Minimum fraction of spots expressing a gene to include it in TESSERA model fitting. |
| `tessera_min_nonz_mean_counts` | `float` | `1.1` | Minimum average count in non-zero expressing spots. |
| `tessera_fdr_threshold` | `float` | `0.05` | FDR threshold for reporting statistically significant differentially expressed genes. |

### Analysis Configurations (`tessera_analyses`)

Defines specific statistical models, spatial subsets, and hypothesis tests.

```yaml
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

### Explanation
* **`label` (`string`):** Descriptive name for plot titles and summary reports.
* **`filter` (`string`):** Logical expression defining the spot population to test:
  * *Gradient example:* Restricts testing to spots along the halo gradient (`> 0.01` and `< 1`), excluding the soma core itself.
  * *Niche example:* Isolates spots in close proximity to motor neurons (`> 0.25`).
* **`design_formula` (`string`):** R formula specifying the fixed effects:
  * `~ 0 + condition + condition:motor_neuron_density`: Fits condition-specific baseline intercepts and condition-specific interaction slopes along the continuous density gradient.
  * `~ 0 + condition`: Standard two-group comparison across all filtered spots.
* **`contrasts` (`map`):** Linear combinations of fitted model coefficients to test:
  * `slope_ALS_vs_CTRL`: Tests whether the spatial distance-decay slope differs between ALS and CTRL.
  * `slope_ALS` / `slope_CTRL`: Tests if individual slopes are significantly non-zero.
  * `ALS_vs_CTRL`: Tests standard differential expression between ALS and CTRL within the filtered niche.

### Outputs (`05_gene_comparison/{analysis}/`)

| Output Target | Description |
| :--- | :--- |
| `objects/tessera_data.rds` | TESSERA data container containing spatial adjacency structures, design matrices, and count lists. |
| `objects/tessera_fits.rds` | List of fitted TESSERA spatial GLMM models for all tested genes. |
| `tables/performance_summary.csv` | Quality metrics (Moran's I autocorrelation of raw counts vs. fitted residuals) per gene and sample. |
| `tables/de_results_all.csv` | Full Wald test statistics, estimates, standard errors, empirical null parameters, and FDR values for all genes. |
| `tables/de_results_significant.csv` | Filtered table of statistically significant genes passing the FDR threshold (`padj < tessera_fdr_threshold`). |
| `plots/volcano/volcano_{contrast}.png` | Volcano plots displaying Log2FC vs. $-log_{10}(\text{adj.\ } P)$ with the top significant genes labeled. |
| `plots/ma/ma_{contrast}.png` | MA plots showing average expression vs. Log2FC for each evaluated contrast. |
| `plots/moran_qc.png` | Boxplot of Moran's I before and after model fitting to confirm removal of spatial autocorrelation. |
| `plots/gradient_profiles/gradient_profile_{contrast}.png` | *(Gradient analyses)* 10-bin mean $\pm$ SE normalized expression curves across the density gradient for top-ranked candidate genes. |
| `plots/heatmap/heatmap_{contrast}.png` | *(Niche analyses)* Sample-by-compartment balanced Z-score expression heatmap for the top 30 candidate genes ranked by adjusted $P$-value. |
