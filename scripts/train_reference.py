import sys
import os
import torch
import rich.pretty
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from matplotlib import rcParams

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import cell2location
from cell2location.models import RegressionModel
from cell2location.utils.filtering import filter_genes

# -----------------------------------------------------------------------------
# Setup & Logging
# -----------------------------------------------------------------------------
if hasattr(snakemake, "log") and snakemake.log:
    log_file = open(snakemake.log[0], "w")
    sys.stdout = log_file
    sys.stderr = log_file

rcParams['pdf.fonttype'] = 42 
torch.set_float32_matmul_precision('medium')

celltype_column = snakemake.params["celltype_column"]
sample_column = snakemake.params["sample_column"]
max_epochs = int(snakemake.params["max_epochs"])
batch_size = int(snakemake.params["batch_size"])
cell_count_cutoff = int(snakemake.params["cell_count_cutoff"])
cell_percentage_cutoff2 = float(snakemake.params["cell_percentage_cutoff2"])
nonz_mean_cutoff = float(snakemake.params["nonz_mean_cutoff"])

# -----------------------------------------------------------------------------
# 1. Load Preprocessed Data
# -----------------------------------------------------------------------------
print("Loading preprocessed single-cell data...")
adata_ref = sc.read_mtx(snakemake.input["mtx"]).T
adata_ref.var_names = pd.read_csv(snakemake.input["features"], header=None)[0].values
adata_ref.obs_names = pd.read_csv(snakemake.input["barcodes"], header=None)[0].values

metadata = pd.read_csv(snakemake.input["metadata"], index_col=0)
adata_ref.obs = metadata

# -----------------------------------------------------------------------------
# 2. Check Signature Separability
# -----------------------------------------------------------------------------
def check_signature_separability(adata, cell_type_key):
    """Compute pairwise signature correlations and warn for high similarity (> 0.95)."""
    cell_types = adata.obs[cell_type_key].unique()

    signatures = {}
    for ct in cell_types:
        mask = adata.obs[cell_type_key] == ct
        X_ct = adata[mask].X
        signatures[ct] = np.asarray(X_ct.mean(axis=0)).flatten()

    has_high_similarity = False
    for ct1 in cell_types:
        for ct2 in cell_types:
            if ct1 < ct2:
                r, _ = pearsonr(signatures[ct1], signatures[ct2])
                if r > 0.95:
                    has_high_similarity = True
                    print(f"WARNING: High signature similarity between '{ct1}' and '{ct2}' (r={r:.3f} > 0.95).")
    
    if not has_high_similarity:
        print("Signature separability checked: all pairwise correlations <= 0.95.")

check_signature_separability(adata_ref, celltype_column)

# -----------------------------------------------------------------------------
# 3. Filter Mitochondrial and Low-Expression Genes
# -----------------------------------------------------------------------------
print("Filtering mitochondrial and lowly expressed genes...")
n_genes_before = adata_ref.shape[1]
n_cells = adata_ref.shape[0]

genes_to_keep = [g for g in adata_ref.var_names if not g.startswith(('MT-', 'mt-'))]
adata_ref = adata_ref[:, genes_to_keep].copy()

selected = filter_genes(
    adata_ref, 
    cell_count_cutoff=cell_count_cutoff, 
    cell_percentage_cutoff2=cell_percentage_cutoff2, 
    nonz_mean_cutoff=nonz_mean_cutoff
)

os.makedirs(os.path.dirname(snakemake.output["filter_plot"]), exist_ok=True)
plt.savefig(snakemake.output["filter_plot"], bbox_inches='tight')
plt.close('all')

adata_ref = adata_ref[:, selected].copy()
n_genes_after = adata_ref.shape[1]
n_genes_filtered = n_genes_before - n_genes_after

print("=" * 50)
print("GENE FILTERING SUMMARY:")
print(f"  - Genes filtered out: {n_genes_filtered}")
print(f"  - Genes remaining:    {n_genes_after}")
print(f"  - Total cells:        {n_cells}")
print("=" * 50)

# -----------------------------------------------------------------------------
# 4. Initialize and Train Regression Model
# -----------------------------------------------------------------------------
print("Setting up RegressionModel...")
cell2location.models.RegressionModel.setup_anndata(
    adata=adata_ref,
    batch_key=sample_column,
    labels_key=celltype_column
)

mod = RegressionModel(adata_ref)
mod.view_anndata_setup()

print(f"Training regression model on GPU (epochs: {max_epochs}, batch size: {batch_size})...")
num_workers = snakemake.threads
mod.train(
    max_epochs=max_epochs, 
    accelerator="gpu", 
    batch_size=batch_size,
    datasplitter_kwargs={'num_workers': num_workers},
)

# Plot training ELBO convergence
os.makedirs(os.path.dirname(snakemake.output["history_plot"]), exist_ok=True)
elbo_df = mod.history["elbo_train"]
loss_col = elbo_df.columns[0]
elbo_df_subset = elbo_df.iloc[20:] if len(elbo_df) > 20 else elbo_df

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(
    elbo_df_subset.index, 
    elbo_df_subset[loss_col], 
    color="#1f77b4",
    linewidth=2,
    label="Training ELBO"
)
ax.set_title("Reference Model Training Convergence", fontsize=12, fontweight='bold')
ax.set_xlabel("Training Steps", fontsize=10)
ax.set_ylabel("ELBO Loss", fontsize=10)
ax.grid(True, linestyle="--", alpha=0.5)
ax.legend(loc="upper right")
fig.savefig(snakemake.output["history_plot"], bbox_inches='tight')
plt.close(fig)

# -----------------------------------------------------------------------------
# 5. Export Posterior and Save Results
# -----------------------------------------------------------------------------
print("Exporting posterior distributions...")
adata_ref = mod.export_posterior(
    adata_ref, 
    sample_kwargs={'num_samples': 1000, 'batch_size': batch_size, 'accelerator': 'gpu'}
)

# Intercept plt.show to capture individual plots from mod.plot_QC()
original_show = plt.show
show_count = 0

def custom_show():
    global show_count
    show_count += 1
    if show_count == 1:
        plt.savefig(snakemake.output["qc_reconstruction_plot"], bbox_inches='tight')
    elif show_count == 2:
        plt.savefig(snakemake.output["qc_expression_plot"], bbox_inches='tight')
    plt.close('all')

plt.show = custom_show
try:
    mod.plot_QC()
finally:
    plt.show = original_show

# Save artifacts
mod.save(snakemake.output["model_dir"], overwrite=True)
adata_ref.write_h5ad(snakemake.output["h5ad"])
print("Reference signature regression script completed successfully.")