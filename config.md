# Pipeline Configuration Guide (`config.yaml`)

This guide details all configuration options for the Visium HD Motor Neuron Spatial Pipeline. Parameters are organized to mirror the workflow execution and output stages.

---

## 1. Global & Input Settings

| Parameter | Type | Default | Rationale & Guidance |
| :--- | :--- | :--- | :--- |
| `output_dir` | `string` | `"FINAL"` | Directory where all rule outputs, intermediate files, models, and plots are stored. |
| `seurat_rds` | `string` | `None` | Path to single-nucleus/single-cell Seurat RDS object containing raw UMI counts and cell type labels. |
| `gm_table_path` | `string` | `None` | Path to a CSV table mapping sample IDs and barcodes belonging to spinal cord grey matter. |
| `visium_samples` | `map` | `None` | Dictionary of Visium HD samples. Each entry specifies the Space Ranger `path` (to `square_016um/`) and sample `condition` (e.g., `ALS`, `CTRL`). |

---

## 2. Stage 01: Single-Cell Reference Prep & Model Training

Prepares reference signatures and fits the negative binomial regression model (`cell2location.models.RegressionModel`).

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

---

## 3. Stage 02: Spatial Deconvolution (cell2location)

Maps reference cell signatures onto 16 µm Visium HD bins (`cell2location.models.Cell2location`).

| Parameter | Type | Default | Rationale & Guidance |
| :--- | :--- | :--- | :--- |
| `umi_min` | `integer` | `50` | Minimum total UMI count per 16 µm bin. Filters out empty glass and damaged tissue. |
| `gene_min` | `integer` | `25` | Minimum unique genes detected per bin. |
| `max_epochs_spatial` | `integer` | `500` | Training iterations for spatial mapping. Check convergence via training history plots. |
| `batch_size_spatial` | `integer` | `16000` | Mini-batch size for spatial optimization. Reduce if encountering CUDA OOM errors. |
| `N_cells_per_location` | `integer` | `1` | Prior mean for the number of cells expected per bin. For 16 µm Visium HD bins, ~1 cell per bin is expected. |
| `detection_alpha` | `float` | `20` | Regularization hyperparameter for spot-to-spot technical sensitivity differences. |

---

## 4. Stage 03: Motor Neuron Segmentation & Halo Gradients

Identifies motor neuron somas using core thresholding, DBSCAN, geodesic expansion, marker gene confirmation, and spatial halo decay.

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

---

## 5. Stage 04: Cell Type Abundance Comparisons

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
* **Key (e.g., `is_grey_matter`, `motor_neuron_density`):** Unique internal identifier for the comparison, used in output file names (e.g., `stats_is_grey_matter.csv`).
* **`label` (`string`):** Clean descriptive title displayed on the generated stacked bar charts and strip plots.
* **`filter` (`string`):** A valid R logical expression evaluated on the spot metadata columns:
  * `"is_grey_matter"`: Evaluates all bins within the anatomically defined grey matter.
  * `"motor_neuron_density > 0.25"`: Restricts analysis to the immediate perineuronal microenvironment (~140 µm radius around motor neurons).

---

## 6. Stage 05: Spatial Differential Expression (TESSERA)

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
