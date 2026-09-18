# ------------------------------------------------------------------------------
# Setup Logging and Snakemake Parameters
# ------------------------------------------------------------------------------
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")

meta_file       <- snakemake@input[["metadata"]]
abund_file      <- snakemake@input[["abundance"]]
out_barchart    <- snakemake@output[["barchart"]]
out_stripplots  <- snakemake@output[["stripplots_dir"]]
out_comb_table  <- snakemake@output[["combined_table"]]

samples_cfg     <- snakemake@params[["visium_samples"]]
bin_comparisons <- snakemake@params[["bin_comparisons"]]

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(rlang)
})

dir.create(dirname(out_barchart), recursive = TRUE, showWarnings = FALSE)
dir.create(out_stripplots, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_comb_table), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Condition & Sample Metadata Harmonization
# ------------------------------------------------------------------------------
sample_condition_raw <- tibble(
  sample    = names(samples_cfg),
  condition = sapply(samples_cfg, function(x) x[["condition"]])
)

cond_levels <- unique(sample_condition_raw$condition)
if ("CTRL" %in% cond_levels) {
  cond_levels <- c("CTRL", setdiff(cond_levels, "CTRL"))
}

sample_counts <- sample_condition_raw %>%
  count(condition, name = "n_samples")

sample_condition_map <- sample_condition_raw %>%
  left_join(sample_counts, by = "condition") %>%
  mutate(
    condition = factor(condition, levels = cond_levels),
    condition_label = paste0(condition, "\n(n = ", n_samples, ")"),
    condition_label = factor(
      condition_label,
      levels = paste0(cond_levels, "\n(n = ", sample_counts$n_samples[match(cond_levels, sample_counts$condition)], ")")
    )
  )

# ------------------------------------------------------------------------------
# 2. Load & Join Tables
# ------------------------------------------------------------------------------
message("Loading spot metadata and abundance tables...")
meta_df  <- read_csv(meta_file, show_col_types = FALSE)
abund_df <- read_csv(abund_file, show_col_types = FALSE)

colnames(meta_df)[1]  <- "spot_id"
colnames(abund_df)[1] <- "spot_id"

# Join on all shared keys (e.g., spot_id, sample, condition) to avoid .x/.y collisions
join_by <- intersect(colnames(meta_df), colnames(abund_df))
combined_df <- inner_join(meta_df, abund_df, by = join_by)

celltype_cols <- names(abund_df)[sapply(abund_df, is.numeric)]
celltype_cols <- setdiff(celltype_cols, c("spot_id", "sample", "in_tissue", "array_row", "array_col"))

message("Identified ", length(celltype_cols), " cell type columns.")

# ------------------------------------------------------------------------------
# 3. Dynamic Filtering & Sample-Level Proportions
# ------------------------------------------------------------------------------
sample_props_list <- list()

for (comp_id in names(bin_comparisons)) {
  comp_spec   <- bin_comparisons[[comp_id]]
  comp_label  <- if (is.list(comp_spec) && !is.null(comp_spec$label)) comp_spec$label else comp_id
  comp_filter <- if (is.list(comp_spec) && !is.null(comp_spec$filter)) comp_spec$filter else comp_spec

  message("Filtering comparison: '", comp_label, "' [", comp_filter, "]")

  sub_df <- tryCatch({
    combined_df %>% filter(!!rlang::parse_expr(comp_filter))
  }, error = function(e) {
    stop("Error applying filter expression '", comp_filter, "': ", e$message)
  })

  if (nrow(sub_df) == 0) {
    warning("No spots matched filter: '", comp_filter, "'. Skipping.")
    next
  }

  comp_props <- sub_df %>%
    group_by(sample) %>%
    summarise(across(all_of(celltype_cols), sum, na.rm = TRUE), .groups = "drop") %>%
    pivot_longer(cols = all_of(celltype_cols), names_to = "celltype", values_to = "abundance") %>%
    # Retain samples with zero matching spots so all cohort samples are represented
    complete(
      sample   = names(samples_cfg),
      celltype = celltype_cols,
      fill     = list(abundance = 0)
    ) %>%
    group_by(sample) %>%
    mutate(sample_prop = if (sum(abundance) > 0) abundance / sum(abundance) else 0) %>%
    ungroup() %>%
    inner_join(sample_condition_map, by = "sample") %>%
    mutate(
      comparison_id = comp_id,
      comparison    = comp_label
    )

  sample_props_list[[comp_id]] <- comp_props
}

if (length(sample_props_list) == 0) {
  stop("No valid data points survived the configured filters.")
}

sample_props_all <- bind_rows(sample_props_list)

# ------------------------------------------------------------------------------
# 4. Coordinate Layout & Group Headers
# ------------------------------------------------------------------------------
comparisons_vec <- unique(sample_props_all$comparison)
conditions_vec  <- levels(sample_condition_map$condition)

layout_list <- list()
header_list <- list()

curr_start  <- 1.00
bar_width   <- 0.44
within_step <- 0.46
between_gap <- 0.40

for (comp in comparisons_vec) {
  comp_cond_x <- numeric(length(conditions_vec))
  for (j in seq_along(conditions_vec)) {
    comp_cond_x[j] <- curr_start + (j - 1) * within_step
  }
  end_x <- comp_cond_x[length(conditions_vec)]

  layout_list[[comp]] <- tibble(
    comparison = comp,
    condition  = factor(conditions_vec, levels = conditions_vec),
    x_pos      = comp_cond_x
  )

  header_list[[comp]] <- tibble(
    comparison = comp,
    x_center   = (curr_start + end_x) / 2,
    x_start    = curr_start - (bar_width / 2) + 0.02,
    x_end      = end_x + (bar_width / 2) - 0.02
  )

  curr_start <- end_x + bar_width + between_gap
}

layout_df <- bind_rows(layout_list)
header_df <- bind_rows(header_list)

sample_props_all <- sample_props_all %>%
  inner_join(layout_df, by = c("comparison", "condition"))

x_axis_ticks <- layout_df %>%
  inner_join(sample_condition_map %>% distinct(condition, condition_label), by = "condition") %>%
  arrange(x_pos)

x_limits <- c(min(x_axis_ticks$x_pos) - 0.55, max(x_axis_ticks$x_pos) + 0.55)

condition_props_all <- sample_props_all %>%
  group_by(comparison_id, comparison, condition, condition_label, x_pos, celltype) %>%
  summarise(mean_prop = mean(sample_prop), .groups = "drop") %>%
  group_by(comparison, condition_label) %>%
  mutate(mean_prop = mean_prop / sum(mean_prop)) %>%
  ungroup()

# ------------------------------------------------------------------------------
# 5. Stacked Bar Chart (Okabe-Ito Palette)
# ------------------------------------------------------------------------------
message("Generating stacked bar chart...")

distinct_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
other_color     <- "#999999"
n_colors        <- length(distinct_colors)

ranked_celltypes <- condition_props_all %>%
  group_by(celltype) %>%
  summarise(global_mean = mean(mean_prop), .groups = "drop") %>%
  arrange(desc(global_mean))

total_celltypes <- nrow(ranked_celltypes)

if (total_celltypes <= n_colors) {
  top_celltypes <- ranked_celltypes$celltype
  plot_bar_df <- condition_props_all %>%
    rename(celltype_plot = celltype) %>%
    mutate(celltype_plot = factor(celltype_plot, levels = top_celltypes))
  color_map <- setNames(distinct_colors[seq_len(total_celltypes)], top_celltypes)
} else {
  top_celltypes <- ranked_celltypes$celltype[seq_len(n_colors)]
  plot_bar_df <- condition_props_all %>%
    mutate(celltype_plot = ifelse(celltype %in% top_celltypes, celltype, "Other")) %>%
    group_by(comparison, condition_label, x_pos, celltype_plot) %>%
    summarise(mean_prop = sum(mean_prop), .groups = "drop") %>%
    mutate(celltype_plot = factor(celltype_plot, levels = c(top_celltypes, "Other")))
  color_map <- setNames(c(distinct_colors, other_color), c(top_celltypes, "Other"))
}

p_bar <- ggplot(plot_bar_df, aes(x = x_pos, y = mean_prop, fill = celltype_plot)) +
  geom_col(width = bar_width, color = "white", linewidth = 0.3) +
  geom_segment(
    data = header_df,
    aes(x = x_start, xend = x_end, y = 1.03, yend = 1.03),
    linewidth = 0.5, color = "grey35", inherit.aes = FALSE
  ) +
  geom_text(
    data = header_df,
    aes(x = x_center, y = 1.07, label = comparison),
    fontface = "bold", size = 3.6, color = "black", angle = 0, inherit.aes = FALSE
  ) +
  scale_y_continuous(
    labels = percent_format(),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_continuous(
    breaks = x_axis_ticks$x_pos,
    labels = x_axis_ticks$condition_label,
    limits = x_limits,
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(0, 1.12), clip = "off") +
  scale_fill_manual(values = color_map) +
  labs(
    title = "Cell Type Abundance Comparison",
    x     = NULL,
    y     = "Cell Type Proportion",
    fill  = "Cell Type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(size = 13, face = "bold", hjust = 0),
    axis.text.x   = element_text(size = 10, face = "bold", color = "black"),
    axis.text.y   = element_text(size = 10, color = "black"),
    axis.title.y  = element_text(size = 11, face = "bold", margin = margin(r = 8)),
    axis.line     = element_line(linewidth = 0.5, color = "grey40"),
    legend.title  = element_text(size = 10.5, face = "bold"),
    legend.text   = element_text(size = 9.5),
    legend.position = "right",
    plot.margin   = margin(t = 18, r = 20, b = 15, l = 15)
  )

bar_plot_width <- max(6.5, 3.5 + 2.0 * length(comparisons_vec))
ggsave(out_barchart, plot = p_bar, width = bar_plot_width, height = 5.5, dpi = 300)
message("Saved stacked bar chart to: ", out_barchart)

# ------------------------------------------------------------------------------
# 6. Stripplots for All Cell Types
# ------------------------------------------------------------------------------
message("Generating individual stripplots for all cell types...")

default_cond_colors <- c("CTRL" = "#0072B2", "ALS" = "#D55E00")
available_conds     <- levels(sample_condition_map$condition)
cond_palette        <- setNames(
  c("#0072B2", "#D55E00", "#009E73", "#E69F00")[seq_along(available_conds)],
  available_conds
)
for (c_name in names(default_cond_colors)) {
  if (c_name %in% available_conds) cond_palette[c_name] <- default_cond_colors[c_name]
}

mean_summary_df <- sample_props_all %>%
  group_by(celltype, comparison, condition, x_pos) %>%
  summarise(mean_prop = mean(sample_prop), .groups = "drop")

plot_celltype_abundance <- function(target_celltype) {
  sub_points <- sample_props_all %>% filter(celltype == target_celltype)
  sub_means  <- mean_summary_df %>% filter(celltype == target_celltype)

  max_val <- max(c(sub_points$sample_prop, sub_means$mean_prop), na.rm = TRUE)
  if (is.na(max_val) || max_val == 0) max_val <- 0.01
  y_top  <- max_val * 1.25
  y_text <- max_val * 1.18
  y_line <- max_val * 1.12

  ggplot() +
    geom_segment(
      data = header_df,
      aes(x = x_start, xend = x_end, y = y_line, yend = y_line),
      linewidth = 0.5, color = "grey35", inherit.aes = FALSE
    ) +
    geom_text(
      data = header_df,
      aes(x = x_center, y = y_text, label = comparison),
      fontface = "bold", size = 3.6, color = "black", angle = 0, inherit.aes = FALSE
    ) +
    geom_col(
      data = sub_means,
      aes(x = x_pos, y = mean_prop, fill = condition, color = condition),
      width = bar_width,
      alpha = 0.35,
      linewidth = 0.5
    ) +
    geom_point(
      data = sub_points,
      aes(x = x_pos, y = sample_prop, fill = condition),
      position = position_jitter(width = 0.04, height = 0, seed = 42),
      shape = 21,
      size = 3.2,
      color = "black",
      stroke = 0.6
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 0.1),
      breaks = pretty_breaks(n = 5),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_x_continuous(
      breaks = x_axis_ticks$x_pos,
      labels = x_axis_ticks$condition_label,
      limits = x_limits,
      expand = c(0, 0)
    ) +
    coord_cartesian(ylim = c(0, y_top), clip = "off") +
    scale_fill_manual(values = cond_palette) +
    scale_color_manual(values = cond_palette) +
    labs(
      title = paste0("Cell Type Abundance: ", target_celltype),
      x     = NULL,
      y     = "Cell Type Proportion"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title    = element_text(size = 13, face = "bold", hjust = 0),
      axis.text.x   = element_text(size = 10, face = "bold", color = "black"),
      axis.text.y   = element_text(size = 10, color = "black"),
      axis.title.y  = element_text(size = 11, face = "bold", margin = margin(r = 8)),
      axis.line     = element_line(linewidth = 0.5, color = "grey40"),
      legend.position = "none",
      plot.margin   = margin(t = 18, r = 20, b = 15, l = 15)
    )
}

stripplot_width <- max(6.0, 3.0 + 2.0 * length(comparisons_vec))

for (ct in celltype_cols) {
  safe_name <- gsub("[^A-Za-z0-9_]+", "_", ct)
  p_ct <- plot_celltype_abundance(ct)
  out_ct_file <- file.path(out_stripplots, paste0("abundance_", safe_name, ".png"))
  ggsave(out_ct_file, plot = p_ct, width = stripplot_width, height = 5.5, dpi = 300)
}
message("Saved all stripplots to: ", out_stripplots)

# ------------------------------------------------------------------------------
# 7. Summary Statistics Tables
# ------------------------------------------------------------------------------
message("Computing summary statistics across all comparisons...")

summary_stats_all <- sample_props_all %>%
  group_by(comparison_id, comparison, celltype, condition) %>%
  summarise(
    mean_prop = mean(sample_prop),
    sd_prop   = ifelse(n() > 1, sd(sample_prop), NA_real_),
    .groups   = "drop"
  ) %>%
  pivot_wider(
    names_from  = condition,
    values_from = c(mean_prop, sd_prop)
  )

ctrl_name <- cond_levels[1]
case_name <- if (length(cond_levels) > 1) cond_levels[2] else cond_levels[1]

ctrl_mean_col <- paste0("mean_prop_", ctrl_name)
case_mean_col <- paste0("mean_prop_", case_name)
ctrl_sd_col   <- paste0("sd_prop_", ctrl_name)
case_sd_col   <- paste0("sd_prop_", case_name)

ctrl_n <- sample_counts$n_samples[sample_counts$condition == ctrl_name]
case_n <- sample_counts$n_samples[sample_counts$condition == case_name]

format_table <- function(df_subset) {
  ctrl_m <- df_subset[[ctrl_mean_col]]
  case_m <- df_subset[[case_mean_col]]
  ctrl_s <- if (ctrl_sd_col %in% names(df_subset)) df_subset[[ctrl_sd_col]] else rep(NA_real_, nrow(df_subset))
  case_s <- if (case_sd_col %in% names(df_subset)) df_subset[[case_sd_col]] else rep(NA_real_, nrow(df_subset))
  l2fc   <- log2(case_m / ctrl_m)

  tibble(
    `Cell Type`                                   = df_subset$celltype,
    !!paste0(ctrl_name, " Mean (n=", ctrl_n, ")") := percent(ctrl_m, accuracy = 0.1),
    !!paste0(ctrl_name, " SD")                    := ifelse(is.na(ctrl_s), "-", percent(ctrl_s, accuracy = 0.1)),
    !!paste0(case_name, " Mean (n=", case_n, ")") := percent(case_m, accuracy = 0.1),
    !!paste0(case_name, " SD")                    := ifelse(is.na(case_s), "-", percent(case_s, accuracy = 0.1)),
    `Log2FC` = case_when(
      is.na(l2fc)       ~ "-",
      is.infinite(l2fc) ~ ifelse(l2fc > 0, "+Inf", "-Inf"),
      TRUE              ~ sprintf("%+.2f", l2fc)
    )
  ) %>%
    arrange(desc(case_m))
}

tables_dir <- dirname(out_comb_table)
all_formatted_list <- list()

for (comp_id in names(bin_comparisons)) {
  comp_spec  <- bin_comparisons[[comp_id]]
  comp_label <- if (is.list(comp_spec) && !is.null(comp_spec$label)) comp_spec$label else comp_id

  sub_stat <- summary_stats_all %>% filter(comparison_id == comp_id)
  if (nrow(sub_stat) == 0) next

  fmt_tbl <- format_table(sub_stat)
  all_formatted_list[[comp_id]] <- fmt_tbl %>% mutate(Comparison = comp_label, .before = 1)

  out_csv <- file.path(tables_dir, paste0("stats_", comp_id, ".csv"))
  write_csv(fmt_tbl, out_csv)
  message("Saved comparison table: ", out_csv)
}

combined_formatted_tbl <- bind_rows(all_formatted_list)
write_csv(combined_formatted_tbl, out_comb_table)
message("Saved master table to: ", out_comb_table)

sink()
sink(type = "message")