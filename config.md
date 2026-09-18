# Pipeline Configuration & Outputs Guide (`config.yaml`)

This guide details all configuration parameters and generated outputs for the Visium HD Motor Neuron Spatial Pipeline, organized to mirror the workflow execution and output stages.

---

## Global & Input Settings

```yaml
output_dir: "FINAL"
seurat_rds: "/path/to/SpinalCord_SingleNucleus_v2.rds"
gm_table_path: "/path/to/Barcodes_GreyMatter.csv"

visium_samples:
  UMC-CD-x030-s:
    path: "/path/to/spaceranger/outs/binned_outputs/square_016um/"
    condition: "ALS"
  UMC-CD-x033-s:
    path: "/path/to/spaceranger/outs/binned_outputs/square_016um/"
    condition: "CTRL"
```

### Explanation & Fields
* **`output_dir` (`string`):** Base directory where all rule outputs, intermediate files, models, tables, and figures are stored (e.g., `"FINAL"`).
* **`seurat_rds` (`string`):** Path to the single-nucleus or single-cell Seurat RDS object containing raw UMI counts and cell type labels.
* **`gm_table_path` (`string`):** Path to a CSV table containing spot-level anatomical classifications. Must contain `sample_id` and `Barcode` columns mapping spots located within the spinal cord grey matter.
* **`visium_samples` (`map`):** Dictionary specifying each Visium HD sample to process:
  * **Sample Key (e.g., `UMC-CD-x030-s`):** Unique identifier used in sample-level file naming.
  * **`path` (`string`):** Directory path to the Space Ranger binned output containing `square_016um/` files (`filtered_feature_bc_matrix.h5`, `tissue_positions.parquet`, etc.).
  * **`condition` (`string`):** Biological cohort or group label (e.g., `ALS`, `CTRL`) used to construct design matrices and contrasts.

---

## Stage 01: Single-Cell Reference Prep & Model Training

Prepares reference signatures and fits the negative binomial regression model (`cell2location.models.RegressionModel`).

### Configuration Parameters

| Parameter | Type | Default | Rationale & Guidance |
| :--- | :--- | :--- | :--- |
| `celltype_column` | `string` | `"motor_neuron_predicted"` | Column in Seurat metadata containing cell type annotations. |
| `group_column` | `string` \| `null`| `"disease_family"` | Metadata column used for balanced subsampling across groups/cohorts. Set to `null` to disable. |
| `sample_column` | `string` | `"orig.ident"` | Metadata column denoting biological library/batch to model technical batch variation. |
| `max_cells_per_celltype` | `integer` | `2000` | Maximum cells sampled per cell type to prevent abundant types from biasing the model or exhausting GPU RAM. |
| `cell_count_cutoff` | `integer` | `5` | Minimum total counts required across cells for a gene to be retained. |
| `cell_percentage_cutoff2` | `float` | `0.05` | Minimum detection rate (5% of cells) in at least one cluster. Retains cluster-specific markers. |
| `nonz_mean_cutoff` | `float` | `1.2` | Minimum average expression of a gene within non-zero expressing cells. |
| `max_epochs` | `integer` | `700` | Total training epochs for the reference regression model. |
| `batch_size` | `integer` | `16000` | Mini-batch size for GPU training. Reduce to `4096` or `2048` if training on smaller GPUs (<40GB VRAM). |

### Stage Outputs (`01_reference/`)

| Output Target | Format | Description |
| :--- | :--- | :--- |
| `tables/cell_counts_summary.csv` | CSV Table | Cell counts per cell type and disease family retained after stratified subsampling. |
| `reference.h5ad` | AnnData H5AD | Downsampled reference AnnData containing exported posterior expression signatures. |
| `models/reference_model/` | Directory | Saved PyTorch/Pyro weights and parameters for the trained cell2location regression model. |
| `plots/filtering_summary.png` | PNG Plot | Diagnostic plot showing gene inclusion based on cell count, non-zero mean, and cluster percentage. |
| `plots/training_history.png` | PNG Plot | ELBO loss convergence curve across reference training epochs. |
| `plots/training_qc_reconstruction.png` | PNG Plot | QC scatter plot evaluating observed vs. model-reconstructed count distributions. |
| `plots/training_qc_expression.png` | PNG Plot | QC scatter plot comparing inferred cluster expression signatures against empirical averages. |

---

## Stage 02: Spatial Deconvolution (cell2location)

Maps reference cell signatures onto 16 µm Visium HD bins (`cell2location.models.Cell2location`).

### Configuration Parameters

| Parameter | Type | Default | Rationale & Guidance |
| :--- | :--- | :--- | :--- |
| `umi_min` | `integer` | `50` | Minimum total UMI count per 16 µm bin. Filters out empty glass and damaged tissue. |
| `gene_min` | `integer` | `25` | Minimum unique genes detected per bin. |
| `max_epochs_spatial` | `integer` | `500` | Training iterations for spatial mapping. Check convergence via training history plots. |
| `batch_size_spatial` | `integer` | `16000` | Mini-batch size for spatial optimization. Reduce if encountering CUDA OOM errors. |
| `N_cells_per_location` | `integer` | `1` | Prior mean for the number of cells expected per bin. For 16 µm Visium HD bins, ~1 cell per bin is expected. |
| `detection_alpha` | `float` | `20` | Regularization hyperparameter for spot-to-spot technical sensitivity differences. |

### Stage Outputs (`02_deconvolution/`)

| Output Target | Format | Description |
| :--- | :--- | :--- |
| `spatial.h5ad` | AnnData H5AD | Combined spatial AnnData containing cell abundance estimates (q05, q50, q95) across all slides. |
| `models/spatial_model/` | Directory | Saved model weights and parameters for the trained spatial cell2location model. |
| `plots/combined_qc_umi.png` | PNG Plot | Histogram of log10 total UMIs across all bins with the cutoff threshold marked. |
| `plots/combined_qc_genes.png` | PNG Plot | Histogram of log10 detected genes across all bins with the cutoff threshold marked. |
| `plots/spatial_training_history.png` | PNG Plot | ELBO loss convergence curve during spatial model training. |
| `plots/spatial_qc_reconstruction.png` | PNG Plot | QC plot evaluating spatial model reconstruction accuracy. |
| `tables/cell_abundance_q05.csv` | CSV Table | 5th-percentile cell abundance estimates per bin for all reference cell types. |
| `tables/spatial_summary_stats.csv` | CSV Table | Summary distribution statistics (mean, SD, median, purity) across spatial bins. |
| `plots/celltypes/*_spatial_abundance.png` | PNG Directory | Multi-slide spatial abundance overlay maps generated individually for each cell type. |

---

## Stage 03: Motor Neuron Segmentation & Halo Gradients

Identifies motor neuron somas using core thresholding, DBSCAN, geodesic expansion, marker gene confirmation, and spatial halo decay.

### Configuration Parameters

| Parameter | Type | Default | Rationale & Guidance |
| :--- | :--- | :--- | :--- |
| `target_col` | `string` | `"Motor Neurons"` | Factor name from cell2location deconvolution representing the target cell type. |
| `bin_size_um` | `float` | `16.0` | Physical width of one Visium HD bin in micrometers. |
| `core_prob_threshold` | `float` | `0.5` | Minimum deconvolution probability required to consider a bin as part of a neuron core. |
| `min_core_bins` | `integer` | `5` | Minimum contiguous bins passing `core_prob_threshold` required to seed a candidate neuron soma. |
| `relaxed_prob_threshold`| `float` | `0.25` | Lower probability boundary for soma expansion from the core. |
| `max_expansion_um` | `float` | `25.0` | Maximum physical radius (in µm) to expand outward from the core via Dijkstra pathfinding. |
| `marker_genes` | `list` | `["CHAT", "SLC5A7"]` | Canonical marker genes required to validate candidate segmented neurons. |
| `min_marker_counts` | `integer` | `1` | Minimum combined raw UMIs of marker genes required across the candidate soma to confirm it. |
| `require_grey_matter` | `boolean` | `true` | When `true`, discards candidate neurons falling outside grey matter to eliminate edge artifacts. |
| `decay_um` | `float` | `100.0` | Distance decay rate ($\lambda$) used in the exponential halo formula: $\text{density} = \exp(-\text{dist} / \lambda)$. |

### Stage Outputs (`03_segmentation/`)

| Output Target | Format | Description |
| :--- | :--- | :--- |
| `tables/motor_neuron_metadata_per_spot.csv` | CSV Table | Spot-level metadata with segmentation flags (`is_neuron_core`, `is_neuron_expanded`, `neuron_id`, `motor_neuron_density`, `is_grey_matter`). |
| `plots/motor_neurons_{sample}_bins.png` | PNG Plot | Histology image showing all verified segmented motor neuron soma bins per sample. |
| `plots/motor_neurons_{sample}_density.png` | PNG Plot | Continuous spatial gradient map illustrating the exponential decay density field ($0 \to 1$) per sample. |

---

## Stage 04: Cell Type Abundance Comparisons

Defines sub-regions or niches where cell type proportions are aggregated and compared between conditions.

```yaml
bin_comparisons:
  is_grey_matter:
    label: "Grey Matter"
    filter: "is_grey_matter"
  motor_neuron_density:
    label: "MN Microenvironment"
    filter: "motor_neuron_density > 0.25"
```

### Explanation & Fields
* **Comparison Key (e.g., `is_grey_matter`, `motor_neuron_density`):** Unique internal identifier for the comparison, used in output file names (e.g., `stats_is_grey_matter.csv`).
* **`label` (`string`):** Clean descriptive title displayed on the generated stacked bar charts and strip plots.
* **`filter` (`string`):** A valid R logical expression evaluated on the spot metadata columns:
  * `"is_grey_matter"`: Evaluates all bins within the anatomically defined grey matter.
  * `"motor_neuron_density > 0.25"`: Restricts analysis to the immediate perineuronal microenvironment (~140 µm radius around motor neurons).

### Stage Outputs (`04_cell_comparison/`)

| Output Target | Format | Description |
| :--- | :--- | :--- |
| `plots/cell_abundance_stacked_barchart.png` | PNG Plot | Stacked bar chart showing condition-averaged cell type proportions across all comparisons. |
| `plots/abundance_stripplots/abundance_*.png` | PNG Directory | Sample-level strip plots with condition mean overlays generated for each individual cell type. |
| `tables/stats_{comparison}.csv` | CSV Table | Condition-level summary metrics (mean proportion, SD, Log2FC) for a specific comparison. |
| `tables/stats_all_comparisons.csv` | CSV Table | Master summary table collating cell type abundance statistics across all configured comparisons. |

---

## Stage 05: Spatial Differential Expression (TESSERA)

Fits spatial generalized linear mixed models (GLMM) with a Leroux CAR random effect to account for spatial autocorrelation across spots.

### Global Parameters

| Parameter | Type | Default | Rationale & Guidance |
| :--- | :--- | :--- | :--- |
| `tessera_d_thresh` | `float` | `23.0` | Maximum pixel distance defining 1st-order orthogonal and diagonal spatial neighbors on the 16 µm grid. |
| `tessera_min_pct_spots` | `float` | `0.01` | Minimum fraction of spots expressing a gene (1%) to include it in TESSERA model fitting. |
| `tessera_min_nonz_mean_counts` | `float` | `1.1` | Minimum average count in non-zero expressing spots. Removes sparse low-count genes. |
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

### Explanation & Fields
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

### Stage Outputs (`05_gene_comparison/{analysis}/`)

| Output Target | Format | Description |
| :--- | :--- | :--- |
| `objects/tessera_data.rds` | RDS Object | TESSERA data container containing spatial adjacency structures, design matrices, and count lists. |
| `objects/tessera_fits.rds` | RDS Object | List of fitted TESSERA spatial GLMM models for all tested genes. |
| `tables/performance_summary.csv` | CSV Table | Quality metrics (Moran's I autocorrelation of raw counts vs. fitted residuals) per gene and sample. |
| `tables/de_results_all.csv` | CSV Table | Full Wald test statistics, estimates, standard errors, empirical null parameters, and FDR values for all genes. |
| `tables/de_results_significant.csv` | CSV Table | Filtered table of statistically significant genes passing the FDR threshold (`padj < tessera_fdr_threshold`). |
| `plots/volcano/volcano_{contrast}.png` | PNG Plot | Volcano plots displaying Log2FC vs. $-log_{10}(\text{adj.\ } P)$ with the top significant genes labeled. |
| `plots/ma/ma_{contrast}.png` | PNG Plot | MA plots showing average expression vs. Log2FC for each evaluated contrast. |
| `plots/moran_qc.png` | PNG Plot | Boxplot of Moran's I before and after model fitting to confirm removal of spatial autocorrelation. |
| `plots/gradient_profiles/gradient_profile_{contrast}.png` | PNG Plot | *(Gradient analyses)* 10-bin mean $\pm$ SE normalized expression curves across the density gradient for top hits. |
| `plots/heatmap/heatmap_{contrast}.png` | PNG Plot | *(Niche analyses)* Sample-by-compartment balanced Z-score expression heatmap for top significant genes. |
