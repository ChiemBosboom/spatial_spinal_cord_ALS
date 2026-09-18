configfile: "config.yaml"

# ==============================================================================
# PATH SHORTCUTS & CONFIG VARIABLES
# ==============================================================================
OUT          = config["output_dir"]
INTERMEDIATE = f"{OUT}/00_intermediate"
REF_DIR      = f"{OUT}/01_reference"
SPATIAL_DIR  = f"{OUT}/02_deconvolution"
SEG_DIR      = f"{OUT}/03_segmentation"
COMP_DIR     = f"{OUT}/04_cell_comparison"
TESSERA_DIR  = f"{OUT}/05_gene_comparison"
BENCHMARKS   = f"{OUT}/benchmarks"
LOGS         = f"{OUT}/logs"

# Wildcard lists
SAMPLES      = list(config["visium_samples"].keys())
COMPARISONS  = list(config["bin_comparisons"].keys())
ANALYSES     = list(config["tessera_analyses"].keys())


# ==============================================================================
# PIPELINE RULES
# ==============================================================================
rule all:
    input:
        # Terminal targets from cell comparison
        f"{COMP_DIR}/tables/stats_all_comparisons.csv",
        f"{COMP_DIR}/plots/cell_abundance_stacked_barchart.png",
        
        # Terminal targets from gene comparison across all analyses
        expand(f"{TESSERA_DIR}/{{analysis}}/tables/de_results_all.csv", analysis=ANALYSES),
        expand(f"{TESSERA_DIR}/{{analysis}}/plots", analysis=ANALYSES)

rule install_tessera:
    output:
        touch(f"{OUT}/.tessera_installed")
    conda:
        "envs/R.yaml"
    log:
        f"{LOGS}/install_tessera.log"
    shell:
        """
        Rscript -e '
        options(repos = c(CRAN = "https://cloud.r-project.org"));
        if (!requireNamespace("sparseinv", quietly = TRUE)) install.packages("sparseinv");
        if (!requireNamespace("TESSERA", quietly = TRUE)) remotes::install_github("floricaconstantine/TESSERA", upgrade = "never", dependencies = FALSE);
        ' > {log} 2>&1
        """

rule downsample_seurat:
    input:
        tessera_ready = f"{OUT}/.tessera_installed",
        seurat_rds = config["seurat_rds"]
    output:
        mtx = temp(f"{INTERMEDIATE}/temp_counts.mtx"),
        features = temp(f"{INTERMEDIATE}/temp_features.tsv"),
        barcodes = temp(f"{INTERMEDIATE}/temp_barcodes.tsv"),
        metadata = temp(f"{INTERMEDIATE}/temp_metadata.csv"),
        summary = f"{REF_DIR}/tables/cell_counts_summary.csv"
    params:
        celltype_column = config["celltype_column"],
        group_column = config["group_column"],
        sample_column = config["sample_column"],
        max_cells_per_celltype = config["max_cells_per_celltype"]
    conda:
        "envs/R.yaml"
    benchmark:
        f"{BENCHMARKS}/downsample_seurat.tsv"
    resources:
        mem_mb = 64000,
        time = "01:00:00"
    threads: 1
    log:
        f"{LOGS}/downsample_seurat.log"
    script:
        "scripts/downsample_seurat.R"

rule train_reference:
    input:
        mtx = rules.downsample_seurat.output.mtx,
        features = rules.downsample_seurat.output.features,
        barcodes = rules.downsample_seurat.output.barcodes,
        metadata = rules.downsample_seurat.output.metadata
    output:
        h5ad = f"{REF_DIR}/reference.h5ad",
        model_dir = directory(f"{REF_DIR}/models/reference_model"),
        history_plot = f"{REF_DIR}/plots/training_history.png",
        qc_reconstruction_plot = f"{REF_DIR}/plots/training_qc_reconstruction.png",
        qc_expression_plot = f"{REF_DIR}/plots/training_qc_expression.png",
        filter_plot = f"{REF_DIR}/plots/filtering_summary.png"
    params:
        celltype_column = config["celltype_column"],
        sample_column = config["sample_column"],
        max_epochs = config["max_epochs"],
        batch_size = config["batch_size"],
        cell_count_cutoff = config["cell_count_cutoff"],
        cell_percentage_cutoff2 = config["cell_percentage_cutoff2"],
        nonz_mean_cutoff = config["nonz_mean_cutoff"]
    conda:
        "envs/python.yaml"
    benchmark:
        f"{BENCHMARKS}/train_reference.tsv"
    resources:
        mem_mb = 128000,
        time = "03:00:00",
        slurm_partition = "gpu",
        gres = "gpu:7g.79gb:1"    
    threads: 4
    log:
        f"{LOGS}/train_reference.log"
    script:
        "scripts/train_reference.py"

rule train_spatial:
    input:
        ref_h5ad = rules.train_reference.output.h5ad,
        ref_model_dir = rules.train_reference.output.model_dir
    output:
        spatial_h5ad = f"{SPATIAL_DIR}/spatial.h5ad",
        spatial_model_dir = directory(f"{SPATIAL_DIR}/models/spatial_model"),
        qc_umi_plot = f"{SPATIAL_DIR}/plots/combined_qc_umi.png",
        qc_genes_plot = f"{SPATIAL_DIR}/plots/combined_qc_genes.png",
        history_plot = f"{SPATIAL_DIR}/plots/spatial_training_history.png",
        spatial_qc_plot = f"{SPATIAL_DIR}/plots/spatial_qc_reconstruction.png",
        abundance_table = f"{SPATIAL_DIR}/tables/cell_abundance_q05.csv",
        summary_table = f"{SPATIAL_DIR}/tables/spatial_summary_stats.csv",
        abundance_plots_dir = directory(f"{SPATIAL_DIR}/plots/celltypes")
    params:
        samples = config["visium_samples"],
        umi_min = config["umi_min"],
        gene_min = config["gene_min"],
        max_epochs_spatial = config["max_epochs_spatial"],
        batch_size_spatial = config["batch_size_spatial"],
        N_cells_per_location = config["N_cells_per_location"],
        detection_alpha = config["detection_alpha"]
    conda:
        "envs/python.yaml"
    benchmark:
        f"{BENCHMARKS}/train_spatial.tsv"
    resources:
        mem_mb = 128000,
        time = "05:00:00",
        slurm_partition = "gpu",
        gres = "gpu:7g.79gb:1"
    threads: 4
    log:
        f"{LOGS}/train_spatial.log"
    script:
        "scripts/train_spatial.py"

rule segment_neurons:
    input:
        spatial_h5ad = rules.train_spatial.output.spatial_h5ad,
        gm_table = config["gm_table_path"]
    output:
        metadata = f"{SEG_DIR}/tables/motor_neuron_metadata_per_spot.csv",
        bins_plots = expand(f"{SEG_DIR}/plots/motor_neurons_{{sample}}_bins.png", sample=SAMPLES),
        density_plots = expand(f"{SEG_DIR}/plots/motor_neurons_{{sample}}_density.png", sample=SAMPLES)
    params:
        target_col = config["target_col"],
        bin_size_um = config["bin_size_um"],
        core_prob_threshold = config["core_prob_threshold"],
        min_core_bins = config["min_core_bins"],
        relaxed_prob_threshold = config["relaxed_prob_threshold"],
        max_expansion_um = config["max_expansion_um"],
        marker_genes = config["marker_genes"],
        min_marker_counts = config["min_marker_counts"],
        require_grey_matter = config["require_grey_matter"],
        decay_um = config["decay_um"]
    conda:
        "envs/python.yaml"
    benchmark:
        f"{BENCHMARKS}/segment_neurons.tsv"
    resources:
        mem_mb = 64000,
        time = "00:45:00"
    threads: 2
    log:
        f"{LOGS}/segment_neurons.log"
    script:
        "scripts/segment_neurons.py"

rule compare_cells:
    input:
        metadata = rules.segment_neurons.output.metadata,
        abundance = rules.train_spatial.output.abundance_table
    output:
        barchart = f"{COMP_DIR}/plots/cell_abundance_stacked_barchart.png",
        stripplots_dir = directory(f"{COMP_DIR}/plots/abundance_stripplots"),
        summary_tables = expand(f"{COMP_DIR}/tables/stats_{{comp}}.csv", comp=COMPARISONS),
        combined_table = f"{COMP_DIR}/tables/stats_all_comparisons.csv"
    params:
        visium_samples = config["visium_samples"],
        bin_comparisons = config["bin_comparisons"]
    conda:
        "envs/R.yaml" 
    benchmark:
        f"{BENCHMARKS}/compare_cells.tsv"
    resources:
        mem_mb = 32000,
        time = "00:30:00"
    threads: 2
    log:
        f"{LOGS}/compare_cells.log"
    script:
        "scripts/compare_cells.R"

rule fit_tessera:
    input:
        metadata = rules.segment_neurons.output.metadata
    output:
        tessera_data = f"{TESSERA_DIR}/{{analysis}}/objects/tessera_data.rds",
        tessera_fits = f"{TESSERA_DIR}/{{analysis}}/objects/tessera_fits.rds",
        perf_summary = f"{TESSERA_DIR}/{{analysis}}/tables/performance_summary.csv"
    params:
        visium_samples = config["visium_samples"],
        umi_min = config["umi_min"],
        gene_min = config["gene_min"],
        d_thresh = config.get("tessera_d_thresh"),
        min_pct_spots = config.get("tessera_min_pct_spots"),
        min_nonz_mean_counts = config.get("tessera_min_nonz_mean_counts"),
        filter_expr = lambda wc: config["tessera_analyses"][wc.analysis]["filter"],
        design_formula = lambda wc: config["tessera_analyses"][wc.analysis]["design_formula"]
    conda:
        "envs/R.yaml"
    benchmark:
        f"{BENCHMARKS}/fit_tessera_{{analysis}}.tsv"
    resources:
        mem_mb = 128000,
        time = "15:00:00"
    threads: 32
    log:
        f"{LOGS}/fit_tessera_{{analysis}}.log"
    script:
        "scripts/fit_tessera.R"

rule compare_genes:
    input:
        tessera_fits = rules.fit_tessera.output.tessera_fits,
        tessera_data = rules.fit_tessera.output.tessera_data
    output:
        results_all = f"{TESSERA_DIR}/{{analysis}}/tables/de_results_all.csv",
        results_sig = f"{TESSERA_DIR}/{{analysis}}/tables/de_results_significant.csv",
        plots_dir = directory(f"{TESSERA_DIR}/{{analysis}}/plots")
    params:
        contrasts = lambda wc: config["tessera_analyses"][wc.analysis].get("contrasts", {}),
        fdr_threshold = config.get("tessera_fdr_threshold", 0.05)
    conda:
        "envs/R.yaml"
    benchmark:
        f"{BENCHMARKS}/compare_genes_{{analysis}}.tsv"
    resources:
        mem_mb = 16000,
        time = "00:30:00"
    threads: 2
    log:
        f"{LOGS}/compare_genes_{{analysis}}.log"
    script:
        "scripts/compare_genes.R"