# ============================================================
# Alpha diversity decomposition by mosquito host
# Effective Shannon, observed richness, and Pielou's evenness
# Total, endogenous, and free viral diversity
# ============================================================

# ----------------------------
# USER SETTINGS
# ----------------------------

alpha_total_file    <- "diversity_results/alpha_total.tsv"
alpha_endo_file     <- "diversity_results/alpha_endo.tsv"
alpha_nonendo_file  <- "diversity_results/alpha_nonendo.tsv"
metadata_file       <- "sample_to_host_clean.csv"

# Use a local output folder if OneDrive causes file-locking issues.
# Example:
# outdir <- "C:/Users/aurowell/Desktop/alpha_diversity_decomposition"

outdir <- "alpha_diversity_decomposition"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# PACKAGES
# ----------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(rstatix)
  library(RColorBrewer)
  library(grid)
})

# ----------------------------
# SAFE WRITE FUNCTION
# ----------------------------

safe_write_csv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  
  tryCatch(
    {
      readr::write_csv(x, path)
    },
    error = function(e) {
      stop(
        "Could not write file:\n",
        path,
        "\n\nLikely causes:\n",
        "- The file is open in Excel\n",
        "- OneDrive is locking or syncing the folder\n",
        "- You do not have write permission to this folder\n\n",
        "Try closing Excel, pausing OneDrive sync, or changing outdir to Desktop/Documents.\n\n",
        "Original error: ",
        conditionMessage(e)
      )
    }
  )
}

# ----------------------------
# LOAD DATA
# ----------------------------

alpha_total <- readr::read_tsv(alpha_total_file, show_col_types = FALSE) |>
  dplyr::mutate(dataset = "Total")

alpha_endo <- readr::read_tsv(alpha_endo_file, show_col_types = FALSE) |>
  dplyr::mutate(dataset = "Endogenous")

alpha_nonendo <- readr::read_tsv(alpha_nonendo_file, show_col_types = FALSE) |>
  dplyr::mutate(dataset = "Free")

meta <- readr::read_csv(metadata_file, show_col_types = FALSE) |>
  dplyr::rename(
    sample = Sample,
    host = Host
  ) |>
  dplyr::distinct(sample, host)

# ----------------------------
# COMBINE DATA
# ----------------------------

alpha_all <- dplyr::bind_rows(alpha_total, alpha_endo, alpha_nonendo) |>
  dplyr::left_join(meta, by = "sample")

if (any(is.na(alpha_all$host))) {
  missing_samples <- alpha_all |>
    dplyr::filter(is.na(host)) |>
    dplyr::distinct(sample) |>
    dplyr::pull(sample)
  
  stop(
    "Some alpha-diversity samples were not found in metadata:\n",
    paste(missing_samples, collapse = ", ")
  )
}

required_cols <- c("sample", "richness", "shannon", "simpson", "dataset", "host")
missing_cols <- setdiff(required_cols, colnames(alpha_all))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns from alpha diversity files:\n",
    paste(missing_cols, collapse = ", ")
  )
}

# ----------------------------
# DERIVED METRICS
# ----------------------------

alpha_all <- alpha_all |>
  dplyr::mutate(
    host = factor(
      host,
      levels = c(
        "Aedes albopictus",
        "Culex pipiens",
        "Culiseta melanura"
      )
    ),
    dataset = factor(
      dataset,
      levels = c("Total", "Endogenous", "Free")
    ),
    effective_shannon = exp(shannon),
    pielou_evenness = dplyr::if_else(
      richness > 1,
      shannon / log(richness),
      NA_real_
    )
  )

# ----------------------------
# LONG FORMAT FOR PLOTTING
# ----------------------------

plot_df <- alpha_all |>
  dplyr::select(
    sample,
    host,
    dataset,
    effective_shannon,
    richness,
    pielou_evenness
  ) |>
  tidyr::pivot_longer(
    cols = c(effective_shannon, richness, pielou_evenness),
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::mutate(
    metric = dplyr::recode(
      metric,
      effective_shannon = "Effective Shannon diversity",
      richness = "Observed richness",
      pielou_evenness = "Pielou's evenness"
    ),
    metric = factor(
      metric,
      levels = c(
        "Effective Shannon diversity",
        "Observed richness",
        "Pielou's evenness"
      )
    )
  )

# ----------------------------
# GLOBAL KRUSKAL-WALLIS TESTS
# ----------------------------

kw_results <- plot_df |>
  dplyr::filter(!is.na(value)) |>
  dplyr::group_by(dataset, metric) |>
  rstatix::kruskal_test(value ~ host) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    p_label = paste0("KW p = ", signif(p, 3))
  )

safe_write_csv(
  kw_results,
  file.path(outdir, "alpha_decomposition_kruskal_results.csv")
)

# ----------------------------
# PAIRWISE WILCOXON TESTS
# ----------------------------

pairwise_results <- plot_df |>
  dplyr::filter(!is.na(value)) |>
  dplyr::group_by(dataset, metric) |>
  rstatix::pairwise_wilcox_test(
    value ~ host,
    p.adjust.method = "BH"
  ) |>
  dplyr::ungroup() |>
  rstatix::add_significance("p.adj")

safe_write_csv(
  pairwise_results,
  file.path(outdir, "alpha_decomposition_pairwise_wilcox_BH.csv")
)

# ----------------------------
# SUMMARY STATISTICS
# ----------------------------

summary_stats <- plot_df |>
  dplyr::group_by(dataset, metric, host) |>
  dplyr::summarise(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    iqr = IQR(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

safe_write_csv(
  summary_stats,
  file.path(outdir, "alpha_decomposition_summary_by_host.csv")
)

# ----------------------------
# PANEL-SPECIFIC Y POSITIONS
# ----------------------------

panel_ranges <- plot_df |>
  dplyr::filter(!is.na(value)) |>
  dplyr::group_by(dataset, metric) |>
  dplyr::summarise(
    ymax = max(value, na.rm = TRUE),
    ymin = min(value, na.rm = TRUE),
    yrange = ymax - ymin,
    .groups = "drop"
  ) |>
  dplyr::mutate(
    yrange = ifelse(yrange == 0, 1, yrange)
  )

# Pairwise brackets:
# - Pielou's evenness brackets are kept within 0-1.
# - Endogenous observed richness is manually lifted so the brackets are legible.
# - Other panels use automatic placement.
pairwise_plot <- pairwise_results |>
  dplyr::left_join(
    panel_ranges,
    by = c("dataset", "metric")
  ) |>
  dplyr::group_by(dataset, metric) |>
  dplyr::mutate(
    y.position = dplyr::case_when(
      metric == "Pielou's evenness" ~ 0.56 + 0.08 * dplyr::row_number(),
      metric == "Observed richness" & dataset == "Endogenous" ~ c(70, 78, 86)[dplyr::row_number()],
      TRUE ~ ymax + yrange * (0.07 + 0.10 * dplyr::row_number())
    )
  ) |>
  dplyr::ungroup()

# KW labels pinned to upper-right of each panel
kw_plot <- kw_results |>
  dplyr::left_join(
    panel_ranges,
    by = c("dataset", "metric")
  ) |>
  dplyr::mutate(
    kw_x = Inf,
    kw_y = Inf
  )

# ----------------------------
# AXIS LIMIT CONTROL
# ----------------------------
# Because facet_grid gives one y-axis scale per row, this section forces
# Pielou's evenness to run from 0 to 1 and gives the other rows enough room
# for brackets and KW labels.

axis_limits_df <- pairwise_plot |>
  dplyr::group_by(metric) |>
  dplyr::summarise(
    max_bracket = max(y.position, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::left_join(
    plot_df |>
      dplyr::group_by(metric) |>
      dplyr::summarise(
        ymin = min(value, na.rm = TRUE),
        ymax = max(value, na.rm = TRUE),
        yrange = ymax - ymin,
        .groups = "drop"
      ),
    by = "metric"
  ) |>
  dplyr::mutate(
    yrange = ifelse(yrange == 0, 1, yrange),
    y_min = dplyr::if_else(metric == "Pielou's evenness", 0, ymin),
    y_max = dplyr::case_when(
      metric == "Pielou's evenness" ~ 1,
      TRUE ~ max_bracket + yrange * 0.08
    ),
    dataset = factor("Total", levels = levels(plot_df$dataset)),
    host = factor("Aedes albopictus", levels = levels(plot_df$host))
  )

axis_limits_long <- axis_limits_df |>
  dplyr::select(dataset, metric, host, y_min, y_max) |>
  tidyr::pivot_longer(
    cols = c(y_min, y_max),
    names_to = "limit_type",
    values_to = "value"
  )

# ----------------------------
# PLOT
# ----------------------------

p <- ggplot(
  plot_df,
  aes(x = host, y = value, color = host)
) +
  geom_blank(
    data = axis_limits_long,
    aes(x = host, y = value),
    inherit.aes = FALSE
  ) +
  geom_boxplot(
    width = 0.55,
    alpha = 0.25,
    outlier.shape = NA,
    linewidth = 0.7
  ) +
  geom_jitter(
    width = 0.12,
    size = 2.3,
    alpha = 0.85
  ) +
  facet_grid(
    rows = vars(metric),
    cols = vars(dataset),
    scales = "free_y",
    switch = "y"
  ) +
  ggpubr::stat_pvalue_manual(
    pairwise_plot,
    label = "p.adj.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.35,
    size = 3.2,
    hide.ns = FALSE
  ) +
  geom_text(
    data = kw_plot,
    aes(
      x = kw_x,
      y = kw_y,
      label = p_label
    ),
    inherit.aes = FALSE,
    hjust = 1.08,
    vjust = 1.35,
    size = 3.2,
    fontface = "bold",
    color = "black"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.03, 0.00))
  ) +
  labs(
    x = "Mosquito host species",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    
    # Column titles
    strip.background.x = element_rect(fill = "gray90", color = NA),
    strip.text.x = element_text(
      face = "bold",
      size = 12,
      margin = margin(t = 4, b = 4)
    ),
    
    # Row titles on y-side
    strip.background.y = element_rect(fill = "gray90", color = NA),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      face = "bold",
      size = 11,
      angle = 90,
      margin = margin(r = 4, l = 4)
    ),
    
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 13),
    
    # Row and panel separators
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.45),
    panel.spacing.x = unit(1.0, "lines"),
    panel.spacing.y = unit(1.1, "lines"),
    
    plot.margin = margin(t = 10, r = 20, b = 15, l = 15)
  )

# ----------------------------
# SAVE OUTPUTS
# ----------------------------

ggsave(
  file.path(outdir, "alpha_decomposition_by_host_3categories_grid.pdf"),
  p,
  width = 16,
  height = 9.5,
  bg = "white"
)

ggsave(
  file.path(outdir, "alpha_decomposition_by_host_3categories_grid.png"),
  p,
  width = 16,
  height = 9.5,
  dpi = 300,
  bg = "white"
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    file.path(outdir, "alpha_decomposition_by_host_3categories_grid.svg"),
    p,
    width = 16,
    height = 9.5,
    device = svglite::svglite,
    bg = "white"
  )
} else {
  warning("Package 'svglite' is not available. SVG file was not saved.")
}

cat("\nAlpha diversity decomposition complete.\n")
cat("Outputs written to:", outdir, "\n")
cat("Main PNG figure:", file.path(outdir, "alpha_decomposition_by_host_3categories_grid.png"), "\n")