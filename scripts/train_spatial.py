import sys
import os
import json
import anndata
import torch
import rich.pretty
import scvi
import pyro
import seaborn as sns
import scanpy as sc
import numpy as np
import pandas as pd
from matplotlib import rcParams
from matplotlib.pyplot import imread
from scvi.module.base import PyroBaseModuleClass

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg')

import cell2location
from cell2location.models import RegressionModel
from cell2location.utils import select_slide

# -----------------------------------------------------------------------------
# Setup & Logging
# -----------------------------------------------------------------------------
if hasattr(snakemake, "log") and snakemake.log:
    log_file = open(snakemake.log[0], "w")
    sys.stdout = log_file
    sys.stderr = log_file

rcParams['pdf.fonttype'] = 42 
torch.set_float32_matmul_precision('medium') 

# -----------------------------------------------------------------------------
# 1. scvi-tools Warmup OOM Patch
# -----------------------------------------------------------------------------
def patched_on_load(self, model, **kwargs):
    """Prevent OOM during scvi-tools load by forcing batch_size=1000 in warmup."""
    pyro.clear_param_store()
    old_history = model.history_.copy() if model.history_ is not None else None
    
    safe_kwargs = self.on_load_kwargs.copy() if self.on_load_kwargs else {}
    safe_kwargs['batch_size'] = 1000 
    
    model.train(max_steps=1, **safe_kwargs)
    model.history_ = old_history
    if "pyro_param_store" in kwargs:
        pyro.get_param_store().set_state(kwargs["pyro_param_store"])

PyroBaseModuleClass.on_load = patched_on_load

# -----------------------------------------------------------------------------
# 2. Visium HD Data Loader
# -----------------------------------------------------------------------------
def read_visium_hd_manual(path, sample_name, condition=None):
    """Load Visium HD 16um bin data with spatial parquet coordinates and scalefactors."""
    print(f"Loading {sample_name} ({condition})...")
    if not path.endswith("/"):
        path = path + "/"

    adata = sc.read_10x_h5(path + "filtered_feature_bc_matrix.h5")
    adata.var_names_make_unique()
    
    positions = pd.read_parquet(path + "spatial/tissue_positions.parquet")
    positions.set_index('barcode', inplace=True)
    
    adata.obs = adata.obs.join(positions)
    if condition is not None:
        adata.obs['condition'] = str(condition)
    
    adata.obsm['spatial'] = adata.obs[['pxl_col_in_fullres', 'pxl_row_in_fullres']].values
    
    with open(path + "spatial/scalefactors_json.json") as f:
        scalefactors = json.load(f)
        
    adata.uns['spatial'] = {
        sample_name: {
            'images': {
                'hires': imread(path + "spatial/tissue_hires_image.png"),
                'lowres': imread(path + "spatial/tissue_lowres_image.png")
            },
            'scalefactors': scalefactors,
            'metadata': {'chemistry_description': "Visium HD", 'software_version': "SpaceRanger"}
        }
    }

    adata.var['SYMBOL'] = adata.var_names
    if 'gene_ids' in adata.var.columns:
        adata.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)
        adata.var_names = adata.var['ENSEMBL']

    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.obs_names = [f"{sample_name}_{i}" for i in adata.obs_names]
    adata.obs.index.name = 'spot_id'

    return adata

# -----------------------------------------------------------------------------
# 3. Concatenate Spatial Slides
# -----------------------------------------------------------------------------
samples = snakemake.params["samples"]
slides = [read_visium_hd_manual(info["path"], name, info.get("condition")) for name, info in samples.items()]

adata_vis = anndata.concat(
    slides, 
    label="sample", 
    keys=list(samples.keys()), 
    uns_merge="unique", 
    merge="same", 
    index_unique=None
)
print(f"Loaded and merged {len(slides)} slides. Total bins: {adata_vis.n_obs}")

# -----------------------------------------------------------------------------
# 4. Generate QC Plots
# -----------------------------------------------------------------------------
UMI_MIN = int(snakemake.params["umi_min"])
GENE_MIN = int(snakemake.params["gene_min"])
plt.style.use('default')

os.makedirs(os.path.dirname(snakemake.output["qc_umi_plot"]), exist_ok=True)

# Plot 1: Combined UMIs
fig1, ax1 = plt.subplots(figsize=(8, 5))
sns.histplot(adata_vis.obs['total_counts'] + 1, log_scale=True, element="step", fill=True, color='blue', ax=ax1)
ax1.axvline(UMI_MIN, color='red', linestyle='--', linewidth=2)
ax1.set_title('Combined Total Counts (UMIs) Across All Samples')
ax1.set_xlabel('UMIs (Log10 Scale)')
ax1.set_ylabel('Number of Bins')
ax1.text(UMI_MIN * 1.1, ax1.get_ylim()[1] * 0.8, f'Min: {UMI_MIN}', color='red', fontweight='bold')
fig1.savefig(snakemake.output["qc_umi_plot"], bbox_inches='tight')
plt.close(fig1)

# Plot 2: Combined Unique Genes
fig2, ax2 = plt.subplots(figsize=(8, 5))
sns.histplot(adata_vis.obs['n_genes_by_counts'] + 1, log_scale=True, element="step", fill=True, color='orange', ax=ax2)
ax2.axvline(GENE_MIN, color='red', linestyle='--', linewidth=2)
ax2.set_title('Combined Unique Genes Across All Samples')
ax2.set_xlabel('Number of Genes (Log10 Scale)')
ax2.set_ylabel('Number of Bins')
ax2.text(GENE_MIN * 1.1, ax2.get_ylim()[1] * 0.8, f'Min: {GENE_MIN}', color='red', fontweight='bold')
fig2.savefig(snakemake.output["qc_genes_plot"], bbox_inches='tight')
plt.close(fig2)

# -----------------------------------------------------------------------------
# 5. Apply Filtering Thresholds
# -----------------------------------------------------------------------------
n_bins_before = adata_vis.n_obs
sc.pp.filter_cells(adata_vis, min_counts=UMI_MIN)
n_after_umi = adata_vis.n_obs

sc.pp.filter_cells(adata_vis, min_genes=GENE_MIN)
n_after_genes = adata_vis.n_obs

print("=" * 50)
print("SPATIAL BIN FILTERING SUMMARY:")
print(f"  - Initial bins:          {n_bins_before}")
print(f"  - Filtered by min UMI:   {n_bins_before - n_after_umi} (Cutoff: {UMI_MIN})")
print(f"  - Filtered by min genes: {n_after_umi - n_after_genes} (Cutoff: {GENE_MIN})")
print(f"  - Bins remaining:        {n_after_genes}")
print("=" * 50)

# -----------------------------------------------------------------------------
# 6. Load Reference and Subset to Shared Genes
# -----------------------------------------------------------------------------
print("Loading reference signatures...")
adata_ref = sc.read_h5ad(snakemake.input["ref_h5ad"])
mod_ref = RegressionModel.load(snakemake.input["ref_model_dir"], adata=adata_ref)

if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']

adata_vis.var_names = adata_vis.var['SYMBOL']
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

print(f"Retained {len(intersect)} overlapping genes between spatial and reference.")

# -----------------------------------------------------------------------------
# 7. Train Cell2location Spatial Model
# -----------------------------------------------------------------------------
print("Setting up cell2location spatial model...")
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

max_epochs_spatial    = int(snakemake.params["max_epochs_spatial"])
batch_size_spatial    = int(snakemake.params["batch_size_spatial"])

N_cells_per_location  = int(snakemake.params["N_cells_per_location"])
detection_alpha       = int(snakemake.params["detection_alpha"])

num_particles_spatial = 4
lr_spatial            = 0.002

mod = cell2location.models.Cell2location(
    adata_vis, 
    cell_state_df=inf_aver,
    N_cells_per_location=N_cells_per_location,
    detection_alpha=detection_alpha
)
mod.view_anndata_setup()

scale_elbo_val = adata_vis.n_obs / batch_size_spatial
print("Training spatial model on GPU...")
mod.train(
    max_epochs=max_epochs_spatial,
    batch_size=batch_size_spatial,
    train_size=1,
    lr=lr_spatial,
    accelerator="gpu",
    num_particles=num_particles_spatial,
    scale_elbo=scale_elbo_val,
    datasplitter_kwargs=dict(
        num_workers=snakemake.threads,
        pin_memory=True
    )
)

# -----------------------------------------------------------------------------
# 8. Save Training History Plot
# -----------------------------------------------------------------------------
elbo_df = mod.history["elbo_train"]
loss_col = elbo_df.columns[0]
elbo_df_subset = elbo_df.iloc[20:] if len(elbo_df) > 20 else elbo_df

fig_hist, ax_hist = plt.subplots(figsize=(7, 4))
ax_hist.plot(elbo_df_subset.index, elbo_df_subset[loss_col], color="#1f77b4", linewidth=2, label="Spatial Mapping ELBO")
ax_hist.set_title("Spatial Model Training Convergence", fontsize=12, fontweight='bold')
ax_hist.set_xlabel("Training Steps", fontsize=10)
ax_hist.set_ylabel("ELBO Loss", fontsize=10)
ax_hist.grid(True, linestyle="--", alpha=0.5)
ax_hist.legend(loc="upper right")

os.makedirs(os.path.dirname(snakemake.output["history_plot"]), exist_ok=True)
fig_hist.savefig(snakemake.output["history_plot"], bbox_inches='tight')
plt.close(fig_hist)

# -----------------------------------------------------------------------------
# 9. Export Posterior & Save Reconstruction QC Plot
# -----------------------------------------------------------------------------
print("Exporting posterior distributions...")
adata = mod.export_posterior(
    adata_vis, 
    use_quantiles=True,
    add_to_obsm=["q05", "q50", "q95"],
    sample_kwargs={'batch_size': 1000, 'accelerator': "gpu"}
)

print("Generating spatial QC reconstruction plot...")
mod.plot_QC(summary_name="q50")
os.makedirs(os.path.dirname(snakemake.output["spatial_qc_plot"]), exist_ok=True)
plt.savefig(snakemake.output["spatial_qc_plot"], bbox_inches="tight")
plt.close('all')

# -----------------------------------------------------------------------------
# 10. Export Cell Abundance Table & Summary Statistics
# -----------------------------------------------------------------------------
factor_names = list(adata.uns['mod']['factor_names'])

print("Exporting cell abundance estimates table...")
abundance_data = adata.obsm['q05_cell_abundance_w_sf']
if isinstance(abundance_data, pd.DataFrame):
    abundance_df = abundance_data.copy()
    abundance_df.columns = factor_names
    abundance_df.index = adata.obs_names
else:
    abundance_df = pd.DataFrame(abundance_data, index=adata.obs_names, columns=factor_names)

abundance_df.insert(0, 'sample', adata.obs['sample'])
if 'condition' in adata.obs.columns:
    abundance_df.insert(1, 'condition', adata.obs['condition'])

abundance_table_path = snakemake.output.get(
    "abundance_table",
    os.path.join(os.path.dirname(snakemake.output["spatial_h5ad"]), "cell_abundance_q05.csv")
)
os.makedirs(os.path.dirname(abundance_table_path), exist_ok=True)
abundance_df.to_csv(abundance_table_path)

# Summary statistics
print("Computing spatial bin summary metrics...")
W = np.asarray(adata.obsm["q05_cell_abundance_w_sf"])
total_cells = W.sum(axis=1)
top1_abundance = W.max(axis=1)
purity = np.divide(top1_abundance, total_cells, out=np.zeros_like(top1_abundance, dtype=float), where=total_cells > 0)

data_rows = {
    "total_cells": total_cells,
    "top_cell_abundance": top1_abundance,
    "purity": purity,
}

summary_stats = [
    {
        "Metric": label,
        "Mean": np.mean(val),
        "Std": np.std(val),
        "Median": np.median(val),
        "Q25": np.percentile(val, 25),
        "Q75": np.percentile(val, 75),
        "Min": np.min(val),
        "Max": np.max(val),
    }
    for label, val in data_rows.items()
]

summary_df = pd.DataFrame(summary_stats).set_index("Metric")
print(summary_df.round(3))

summary_table_path = snakemake.output.get(
    "summary_table",
    os.path.join(os.path.dirname(snakemake.output["spatial_h5ad"]), "spatial_summary_stats.csv")
)
os.makedirs(os.path.dirname(summary_table_path), exist_ok=True)
summary_df.to_csv(summary_table_path)

# -----------------------------------------------------------------------------
# 11. Plot Spatial Abundance per Cell Type
# -----------------------------------------------------------------------------
print("Generating multi-slide spatial plots per cell type...")
adata.obs[factor_names] = adata.obsm['q05_cell_abundance_w_sf']

sorted_sample_items = sorted(
    samples.items(), 
    key=lambda item: (str(item[1].get("condition", "")), item[0])
)
sorted_sample_names = [item[0] for item in sorted_sample_items]

n_samples = len(sorted_sample_names)
ncols = min(2, n_samples)
nrows = int(np.ceil(n_samples / ncols))

abundance_plots_dir = snakemake.output.get(
    "abundance_plots_dir",
    os.path.join(os.path.dirname(snakemake.output["history_plot"]), "celltypes")
)
os.makedirs(abundance_plots_dir, exist_ok=True)

for cell_type in factor_names:
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 5.0 * nrows), squeeze=False)
    
    with mpl.rc_context({'axes.facecolor': 'black', 'axes.grid': False, 'figure.facecolor': 'white'}):
        for idx, (s_name, s_info) in enumerate(sorted_sample_items):
            r = idx // ncols
            c = idx % ncols
            ax = axes[r, c]
            
            cond = s_info.get("condition", "N/A")
            slide = select_slide(adata, s_name)
            
            sc.pl.spatial(
                slide,
                library_id=s_name,
                cmap='plasma',
                color=cell_type,
                size=1,
                img_key='hires',
                alpha=0.5,
                bw=True,
                title=f"{s_name} ({cond})\n{cell_type}",
                ax=ax,
                show=False
            )

        for idx in range(n_samples, nrows * ncols):
            r = idx // ncols
            c = idx % ncols
            axes[r, c].set_visible(False)
            
    fig.tight_layout()
    clean_filename = "".join(c if c.isalnum() or c in ('_', '-') else '_' for c in cell_type)
    fig.savefig(os.path.join(abundance_plots_dir, f"{clean_filename}_spatial_abundance.png"), dpi=500, bbox_inches='tight')
    plt.close(fig)

# -----------------------------------------------------------------------------
# 12. Save Trained Model & AnnData Object
# -----------------------------------------------------------------------------
print("Saving model directory and spatial AnnData...")
mod.save(snakemake.output["spatial_model_dir"], overwrite=True)
adata.write_h5ad(snakemake.output["spatial_h5ad"])
print("Spatial training script completed successfully.")