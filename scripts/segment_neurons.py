import os
import sys
import re
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import issparse
from sklearn.cluster import DBSCAN

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Setup & Logging
# -----------------------------------------------------------------------------
if hasattr(snakemake, "log") and snakemake.log:
    log_file = open(snakemake.log[0], "w")
    sys.stdout = log_file
    sys.stderr = log_file

matplotlib.rcParams["pdf.fonttype"] = 42

spatial_h5ad_path   = snakemake.input["spatial_h5ad"]
gm_table_path       = snakemake.input["gm_table"]

output_metadata_csv = snakemake.output["metadata"]
bins_plot_files     = snakemake.output["bins_plots"]
density_plot_files  = snakemake.output["density_plots"]

target_col             = snakemake.params["target_col"]
bin_size_um            = snakemake.params["bin_size_um"]
core_prob_threshold    = snakemake.params["core_prob_threshold"]
min_core_bins          = snakemake.params["min_core_bins"]
relaxed_prob_threshold = snakemake.params["relaxed_prob_threshold"]
max_expansion_um       = snakemake.params["max_expansion_um"]
marker_genes           = snakemake.params["marker_genes"]
min_marker_counts      = snakemake.params["min_marker_counts"]
require_grey_matter    = snakemake.params["require_grey_matter"]
decay_um               = snakemake.params["decay_um"]

# Plotting styling
bin_color    = "#00E5FF"
spot_size    = 1
density_cmap = "viridis"
density_vmin = 0.0
density_vmax = 1.0

def get_sample_filepath(file_list, sample_id):
    """Locate output plot path corresponding to a given sample identifier."""
    for f in file_list:
        if sample_id in f:
            return f
    raise ValueError(f"Could not find output file matching sample '{sample_id}' in {file_list}")

# -----------------------------------------------------------------------------
# 1. Load Data & Map Grey Matter Barcodes
# -----------------------------------------------------------------------------
print(f"Loading spatial data: {spatial_h5ad_path}")
adata_vis = sc.read_h5ad(spatial_h5ad_path)

print(f"Loading Grey Matter table: {gm_table_path}")
gm_df = pd.read_csv(gm_table_path)

extracted_barcodes = adata_vis.obs_names.str.extract(r"(s_\d+um_\d+_\d+-\d+)", expand=False)
valid_samples = gm_df["sample_id"].astype(str).unique()
sample_pattern = r"(" + "|".join(map(re.escape, valid_samples)) + r")"
extracted_samples = adata_vis.obs_names.str.extract(sample_pattern, expand=False)

gm_lookup_set = set(gm_df["sample_id"].astype(str) + "__" + gm_df["Barcode"].astype(str))
obs_composite_keys = extracted_samples + "__" + extracted_barcodes
adata_vis.obs["is_grey_matter"] = obs_composite_keys.isin(gm_lookup_set)

print(f"Mapped {adata_vis.obs['is_grey_matter'].sum():,} grey matter bins.")

# -----------------------------------------------------------------------------
# 2. Extract Marker Expression and Abundance Score
# -----------------------------------------------------------------------------
def extract_gene_expression(adata, gene_list):
    """Aggregate expression for specified marker genes across all bins."""
    combined_expr = np.zeros(adata.n_obs, dtype=np.float32)
    found_genes = []
    source = adata.raw if adata.raw is not None else adata
    var_names_upper = {v.upper(): v for v in source.var_names}

    for gene in gene_list:
        gene_upper = gene.upper()
        if gene_upper in var_names_upper:
            actual_name = var_names_upper[gene_upper]
            found_genes.append(actual_name)
            expr = source[:, actual_name].X
            if issparse(expr):
                expr = expr.toarray()
            combined_expr += expr.ravel()

    print(f"Marker verification initialized using: {found_genes}")
    return combined_expr

adata_vis.obs["MN_Marker_Expr"] = extract_gene_expression(adata_vis, marker_genes)

col_name = f"q05cell_abundance_w_sf_{target_col}" if not target_col.startswith("q05") else target_col
adata_vis.obs["Motor_Neuron_Score"] = adata_vis.obsm["q05_cell_abundance_w_sf"][col_name]

adata_vis.obs["is_neuron_core"] = False
adata_vis.obs["is_neuron_expanded"] = False
adata_vis.obs["neuron_id"] = "None"

px_per_um_dict = {}
unique_samples = adata_vis.obs["sample"].unique()

# -----------------------------------------------------------------------------
# 3. Instance Segmentation & Dual Validation
# -----------------------------------------------------------------------------
for sample in unique_samples:
    sample_mask = adata_vis.obs["sample"] == sample
    sample_obs = adata_vis.obs[sample_mask]
    coords = sample_obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values

    # Scale calibration
    tree_all = cKDTree(coords)
    dists, _ = tree_all.query(coords, k=2)
    px_per_bin = np.median(dists[:, 1])
    px_per_um = px_per_bin / bin_size_um
    px_per_um_dict[sample] = px_per_um

    adj_radius_px = px_per_bin * 1.5
    max_expansion_px = max_expansion_um * px_per_um

    # Core Identification
    core_mask = sample_obs["Motor_Neuron_Score"] >= core_prob_threshold
    if not core_mask.any():
        print(f"  Sample {sample}: No core candidates detected.")
        continue

    core_coords = coords[core_mask]
    core_indices = sample_obs.index[core_mask].values

    db = DBSCAN(eps=adj_radius_px, min_samples=1).fit(core_coords)
    labels, counts = np.unique(db.labels_, return_counts=True)

    valid_labels = labels[counts >= min_core_bins]
    if len(valid_labels) == 0:
        print(f"  Sample {sample}: No core clusters met minimum bin count ({min_core_bins}).")
        continue

    core_is_valid = np.isin(db.labels_, valid_labels)
    valid_core_indices = core_indices[core_is_valid]
    valid_core_coords = core_coords[core_is_valid]
    valid_core_labels = db.labels_[core_is_valid]
    core_neuron_ids = np.array([f"{sample}_Neuron_{lbl}" for lbl in valid_core_labels])

    # Connected Expansion via Dijkstra
    cand_mask = sample_obs["Motor_Neuron_Score"] >= relaxed_prob_threshold
    cand_indices = sample_obs.index[cand_mask].values
    cand_coords = coords[cand_mask]

    tree_cands = cKDTree(cand_coords)
    adj_matrix = tree_cands.sparse_distance_matrix(tree_cands, max_distance=adj_radius_px)

    _, core_in_cand_idx = tree_cands.query(valid_core_coords, k=1)
    shortest_paths = dijkstra(csgraph=adj_matrix, directed=False, indices=core_in_cand_idx)

    min_distances = np.min(shortest_paths, axis=0)
    closest_core_pt = np.argmin(shortest_paths, axis=0)

    reachable = min_distances <= max_expansion_px
    assigned_cand_indices = cand_indices[reachable]
    assigned_neuron_ids = core_neuron_ids[closest_core_pt[reachable]]

    # Marker & Grey Matter Verification
    df_bins = pd.DataFrame({
        "bin_idx": assigned_cand_indices,
        "neuron_id": assigned_neuron_ids,
        "marker_expr": adata_vis.obs.loc[assigned_cand_indices, "MN_Marker_Expr"].values,
        "is_grey_matter": adata_vis.obs.loc[assigned_cand_indices, "is_grey_matter"].values
    })

    marker_sums = df_bins.groupby("neuron_id")["marker_expr"].sum()
    marker_pass = marker_sums >= min_marker_counts

    gm_sums = df_bins.groupby("neuron_id")["is_grey_matter"].sum()
    gm_pass = gm_sums >= 1 if require_grey_matter else pd.Series(True, index=gm_sums.index)

    verified_neuron_mask = marker_pass & gm_pass
    verified_neuron_ids = set(verified_neuron_mask[verified_neuron_mask].index)

    failed_markers = set(marker_pass[~marker_pass].index)
    failed_gm = set(gm_pass[~gm_pass].index)

    if len(verified_neuron_ids) > 0:
        verified_bins_mask = df_bins["neuron_id"].isin(verified_neuron_ids)
        verified_cand_indices = df_bins.loc[verified_bins_mask, "bin_idx"].values
        verified_cand_nids = df_bins.loc[verified_bins_mask, "neuron_id"].values

        adata_vis.obs.loc[verified_cand_indices, "is_neuron_expanded"] = True
        adata_vis.obs.loc[verified_cand_indices, "neuron_id"] = verified_cand_nids

        core_df = pd.DataFrame({"idx": valid_core_indices, "neuron_id": core_neuron_ids})
        verified_core_indices = core_df.loc[core_df["neuron_id"].isin(verified_neuron_ids), "idx"].values
        adata_vis.obs.loc[verified_core_indices, "is_neuron_core"] = True

    print(f"  Sample {sample}: {len(valid_labels)} candidates -> "
          f"{len(verified_neuron_ids)} verified | "
          f"Filtered: {len(failed_markers)} (marker) / {len(failed_gm)} (outside GM)")

# -----------------------------------------------------------------------------
# 4. Spatial Density Gradient Calculation
# -----------------------------------------------------------------------------
print(f"Calculating spatial density gradient (decay = {decay_um} µm)...")
density_col = "motor_neuron_density"
adata_vis.obs[density_col] = 0.0

for sample in unique_samples:
    sample_mask = adata_vis.obs["sample"] == sample
    sample_obs = adata_vis.obs[sample_mask]

    source_mask_sample = sample_obs["is_neuron_expanded"]
    gm_mask_sample = sample_obs["is_grey_matter"]

    if not source_mask_sample.any() or not gm_mask_sample.any():
        print(f"  Sample {sample}: Insufficient neuron sources or GM bins. Skipping density.")
        continue

    coords = sample_obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values
    source_coords = coords[source_mask_sample]
    gm_coords = coords[gm_mask_sample]
    px_per_um = px_per_um_dict[sample]

    tree_sources = cKDTree(source_coords)
    min_dists_px, _ = tree_sources.query(gm_coords, k=1)
    min_dists_um = min_dists_px / px_per_um

    densities_gm = np.exp(-min_dists_um / decay_um)
    target_indices = sample_obs.index[gm_mask_sample]
    adata_vis.obs.loc[target_indices, density_col] = densities_gm

    print(f"  Sample {sample}: Density computed for {len(gm_coords):,} GM bins "
          f"(Mean: {densities_gm.mean():.3f}, Max: {densities_gm.max():.3f}).")

# -----------------------------------------------------------------------------
# 5. Generate Figures
# -----------------------------------------------------------------------------
print("Generating segmentation and density gradient plots...")

for sample in unique_samples:
    sample_mask = adata_vis.obs["sample"] == sample
    sample_obs = adata_vis.obs[sample_mask]

    he_img = None
    scale_factor = 1.0
    if "spatial" in adata_vis.uns and sample in adata_vis.uns["spatial"]:
        spatial_data = adata_vis.uns["spatial"][sample]
        if "hires" in spatial_data.get("images", {}):
            he_img = spatial_data["images"]["hires"]
            scale_factor = spatial_data["scalefactors"]["tissue_hires_scalef"]
        elif "fullres" in spatial_data.get("images", {}):
            he_img = spatial_data["images"]["fullres"]
            scale_factor = 1.0

    # Plot 1: Verified Neuron Bins
    bin_out_path = get_sample_filepath(bins_plot_files, sample)
    os.makedirs(os.path.dirname(bin_out_path), exist_ok=True)

    fig_bins, ax_bins = plt.subplots(figsize=(12, 12), dpi=300)
    if he_img is not None:
        ax_bins.imshow(he_img)
    else:
        ax_bins.set_facecolor("#222222")
        ax_bins.invert_yaxis()

    is_expanded = sample_obs["is_neuron_expanded"].astype(bool)
    is_valid_id = sample_obs["neuron_id"].astype(str) != "None"
    plot_mask_bins = is_expanded & is_valid_id

    if plot_mask_bins.any():
        subset_bins = sample_obs[plot_mask_bins]
        bx = subset_bins["pxl_col_in_fullres"].values * scale_factor
        by = subset_bins["pxl_row_in_fullres"].values * scale_factor
        ax_bins.scatter(
            bx, by,
            c=bin_color,
            s=spot_size,
            edgecolors="black",
            linewidths=0.2,
            zorder=2
        )

    ax_bins.set_title(f"Sample: {sample} | Verified Motor Neuron Bins", fontsize=14, pad=10)
    ax_bins.axis("off")
    fig_bins.tight_layout()
    fig_bins.savefig(bin_out_path, dpi=300, bbox_inches="tight")
    plt.close(fig_bins)

    # Plot 2: Density Gradient
    density_out_path = get_sample_filepath(density_plot_files, sample)
    os.makedirs(os.path.dirname(density_out_path), exist_ok=True)

    fig_dens, ax_dens = plt.subplots(figsize=(12, 12), dpi=300)
    if he_img is not None:
        ax_dens.imshow(he_img)
    else:
        ax_dens.set_facecolor("#222222")
        ax_dens.invert_yaxis()

    is_gm = sample_obs["is_grey_matter"].astype(bool)
    density_vals = sample_obs[density_col].astype(float)
    plot_mask_dens = is_gm & (density_vals > 0.0)

    if plot_mask_dens.any():
        subset_dens = sample_obs[plot_mask_dens]
        dx = subset_dens["pxl_col_in_fullres"].values * scale_factor
        dy = subset_dens["pxl_row_in_fullres"].values * scale_factor
        d_vals = subset_dens[density_col].astype(float).values

        scatter = ax_dens.scatter(
            dx, dy,
            c=d_vals,
            cmap=density_cmap,
            vmin=density_vmin,
            vmax=density_vmax,
            s=spot_size,
            edgecolors="none",
            zorder=2
        )
        cbar = fig_dens.colorbar(scatter, ax=ax_dens, fraction=0.035, pad=0.02, shrink=0.7)
        cbar.set_label("Motor Neuron Density Score", fontsize=11, labelpad=10)
        cbar.ax.tick_params(labelsize=9)

    ax_dens.set_title(f"Sample: {sample} | Motor Neuron Density Gradient ({int(decay_um)}µm decay)", fontsize=14, pad=10)
    ax_dens.axis("off")
    fig_dens.tight_layout()
    fig_dens.savefig(density_out_path, dpi=300, bbox_inches="tight")
    plt.close(fig_dens)

# -----------------------------------------------------------------------------
# 6. Export Spot-Level Metadata Table
# -----------------------------------------------------------------------------
print(f"Exporting metadata table to: {output_metadata_csv}")
os.makedirs(os.path.dirname(output_metadata_csv), exist_ok=True)

desired_columns = [
    "sample",
    "condition",
    "is_grey_matter",
    "is_neuron_core",
    "is_neuron_expanded",
    "neuron_id",
    "motor_neuron_density"
]

export_columns = [col for col in desired_columns if col in adata_vis.obs.columns]
metadata_df = adata_vis.obs[export_columns].copy()
metadata_df.index.name = "spot_id"

for col in ["is_grey_matter", "is_neuron_core", "is_neuron_expanded"]:
    if col in metadata_df.columns:
        metadata_df[col] = metadata_df[col].astype(bool)

if "motor_neuron_density" in metadata_df.columns:
    metadata_df["motor_neuron_density"] = metadata_df["motor_neuron_density"].astype(float)

metadata_df.to_csv(output_metadata_csv, index=True)
print(f"Metadata exported ({metadata_df.shape[0]:,} spots x {metadata_df.shape[1]} columns).")