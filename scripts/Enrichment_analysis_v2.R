# ============================================================
# Enrichment analysis script
# ============================================================

# ------------------------------------------------------------
# Load required packages
# ------------------------------------------------------------

required_packages <- c("tidyverse", "rstatix", "ggpubr", "coin", "openxlsx")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
  library(ggpubr)
})

# Do not load coin with library(coin).
# rstatix::wilcox_effsize() may need coin installed, but loading coin can mask
# rstatix/dplyr functions and cause method conflicts.

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------

input_file <- "summary.csv"

shared_folder <- "enrichment_analysis_results"

output_dir <- shared_folder
figures_dir <- file.path(output_dir, "figures")
tables_dir <- file.path(output_dir, "tables")
processed_dir <- file.path(output_dir, "processed_data")
logs_dir <- file.path(output_dir, "logs")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

safe_log10_plus1 <- function(x) {
  dplyr::if_else(
    is.na(x) | x < 0,
    NA_real_,
    log10(x + 1)
  )
}

safe_percent_change <- function(enriched_value, unenriched_value) {
  dplyr::case_when(
    is.na(unenriched_value) | is.na(enriched_value) ~ NA_real_,
    unenriched_value == 0 & enriched_value == 0 ~ 0,
    unenriched_value == 0 & enriched_value > 0 ~ NA_real_,
    TRUE ~ ((enriched_value - unenriched_value) / unenriched_value) * 100
  )
}

safe_fold_change <- function(enriched_value, unenriched_value) {
  dplyr::case_when(
    is.na(unenriched_value) | is.na(enriched_value) ~ NA_real_,
    unenriched_value == 0 & enriched_value == 0 ~ 1,
    unenriched_value == 0 & enriched_value > 0 ~ NA_real_,
    TRUE ~ enriched_value / unenriched_value
  )
}

format_stat_label <- function(signif, r) {
  paste0(signif, "   r = ", sprintf("%.2f", r))
}

add_stat_positions <- function(plot_df, stats_df, facet_col = "Metric_label") {
  y_positions <- plot_df %>%
    dplyr::group_by(.data[[facet_col]]) %>%
    dplyr::summarise(
      y_min = min(Value, na.rm = TRUE),
      y_max = max(Value, na.rm = TRUE),
      y_range = y_max - y_min,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      y_range = dplyr::case_when(
        is.na(y_range) ~ 1,
        y_range == 0 & !is.na(y_max) ~ abs(y_max) * 0.1,
        TRUE ~ y_range
      ),
      y_range = ifelse(y_range == 0 | is.na(y_range), 1, y_range),
      y_pos = y_max + 0.16 * y_range
    )
  
  stats_df %>%
    dplyr::left_join(y_positions, by = facet_col)
}

# ------------------------------------------------------------
# Read and clean input data
# ------------------------------------------------------------

df <- read.csv(input_file, check.names = FALSE)

names(df) <- names(df) %>%
  stringr::str_replace_all(" ", "_") %>%
  stringr::str_replace_all("\\.", "_") %>%
  stringr::str_replace_all("-", "_")

# ------------------------------------------------------------
# Validate required columns
# ------------------------------------------------------------

required_cols <- c(
  "Sample",
  "Host",
  "Mosquito_reads",
  "Percent_mosquito",
  "Non_mosquito_reads",
  "Total_reads",
  "Viral_reads",
  "Percent_viral",
  "Average_coverage",
  "Mosquito_genome_coverage_per_million_mosquito_reads",
  "Enriched"
)

missing_cols <- setdiff(required_cols, names(df))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nAvailable columns are: ",
      paste(names(df), collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# Prepare variables and derived read metrics
# ------------------------------------------------------------

df_clean <- df %>%
  dplyr::mutate(
    Enriched = dplyr::case_when(
      Enriched %in% c(TRUE, "TRUE", "True", "true", "T", "t", "1", 1) ~ TRUE,
      Enriched %in% c(FALSE, "FALSE", "False", "false", "F", "f", "0", 0) ~ FALSE,
      TRUE ~ NA
    ),
    
    Pair_ID = Sample %>%
      stringr::str_remove("^Sample_") %>%
      stringr::str_remove("^E") %>%
      stringr::str_extract("^[0-9]+"),
    
    Mosquito_reads = as.numeric(Mosquito_reads),
    Viral_reads = as.numeric(Viral_reads),
    Non_mosquito_reads = as.numeric(Non_mosquito_reads),
    Total_reads = as.numeric(Total_reads),
    Percent_mosquito = as.numeric(Percent_mosquito),
    Percent_viral = as.numeric(Percent_viral),
    Average_coverage = as.numeric(Average_coverage),
    Mosquito_genome_coverage_per_million_mosquito_reads =
      as.numeric(Mosquito_genome_coverage_per_million_mosquito_reads),
    
    Normalized_mosquito_coverage =
      Mosquito_genome_coverage_per_million_mosquito_reads,
    
    Percent_mosquito = dplyr::if_else(
      Percent_mosquito > 1,
      Percent_mosquito / 100,
      Percent_mosquito
    ),
    
    Percent_viral = dplyr::if_else(
      Percent_viral > 1,
      Percent_viral / 100,
      Percent_viral
    ),
    
    Other_reads = Non_mosquito_reads - Viral_reads,
    
    Percent_other = Other_reads / Total_reads,
    Percent_non_mosquito = Non_mosquito_reads / Total_reads,
    
    Viral_RPM = Viral_reads / Total_reads * 1e6,
    Other_RPM = Other_reads / Total_reads * 1e6,
    Non_mosquito_RPM = Non_mosquito_reads / Total_reads * 1e6,
    Mosquito_RPM = Mosquito_reads / Total_reads * 1e6,
    
    Viral_percent_non_mosquito = Viral_reads / Non_mosquito_reads,
    Viral_RPM_non_mosquito = Viral_reads / Non_mosquito_reads * 1e6,
    
    log10_Viral_Reads = safe_log10_plus1(Viral_reads),
    log10_Other_Reads = safe_log10_plus1(Other_reads),
    log10_Non_mosquito_Reads = safe_log10_plus1(Non_mosquito_reads),
    log10_Mosquito_Reads = safe_log10_plus1(Mosquito_reads),
    log10_Total_Reads = safe_log10_plus1(Total_reads),
    log10_Viral_RPM = safe_log10_plus1(Viral_RPM),
    log10_Other_RPM = safe_log10_plus1(Other_RPM),
    log10_Viral_RPM_non_mosquito = safe_log10_plus1(Viral_RPM_non_mosquito),
    
    Enrichment = factor(
      Enriched,
      levels = c(FALSE, TRUE),
      labels = c("FALSE", "TRUE")
    ),
    
    Host = factor(
      Host,
      levels = c(
        "Aedes albopictus",
        "Culex pipiens",
        "Culiseta melanura"
      )
    )
  )

if (any(is.na(df_clean$Enriched))) {
  warning("Some Enriched values could not be converted to TRUE/FALSE.")
}

if (any(is.na(df_clean$Pair_ID))) {
  warning("Some Pair_ID values could not be extracted from Sample names.")
}

if (any(is.na(df_clean$Host))) {
  warning("Some Host values did not match the expected host names.")
}

if (any(df_clean$Other_reads < 0, na.rm = TRUE)) {
  warning(
    "Some Other_reads values are negative. Check whether Non_mosquito_reads includes Viral_reads. Negative Other_reads values were set to NA for log10-transformed Other metrics."
  )
}

write.csv(
  df_clean,
  file.path(processed_dir, "cleaned_enrichment_data.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Identify complete enriched/unenriched pairs
# ------------------------------------------------------------

pair_check <- df_clean %>%
  dplyr::filter(!is.na(Pair_ID), !is.na(Enrichment)) %>%
  dplyr::count(Pair_ID, Enrichment, name = "n") %>%
  tidyr::pivot_wider(
    names_from = Enrichment,
    values_from = n,
    values_fill = 0
  )

if (!"FALSE" %in% names(pair_check)) {
  pair_check[["FALSE"]] <- 0
}

if (!"TRUE" %in% names(pair_check)) {
  pair_check[["TRUE"]] <- 0
}

pair_check <- pair_check %>%
  dplyr::mutate(
    complete_pair = `FALSE` == 1 & `TRUE` == 1
  )

write.csv(
  pair_check,
  file.path(tables_dir, "pair_check.csv"),
  row.names = FALSE
)

complete_pairs <- pair_check %>%
  dplyr::filter(complete_pair) %>%
  dplyr::pull(Pair_ID)

if (length(complete_pairs) == 0) {
  stop("No complete enriched/unenriched pairs were detected.")
}

df_paired <- df_clean %>%
  dplyr::filter(Pair_ID %in% complete_pairs) %>%
  dplyr::arrange(Pair_ID, Enrichment)

# ------------------------------------------------------------
# Define metrics for analysis
# ------------------------------------------------------------

metrics_to_analyze <- c(
  "Mosquito_reads",
  "Viral_reads",
  "Other_reads",
  "Non_mosquito_reads",
  "Total_reads",
  "Percent_mosquito",
  "Percent_viral",
  "Percent_other",
  "Percent_non_mosquito",
  "Viral_RPM",
  "Other_RPM",
  "Non_mosquito_RPM",
  "Mosquito_RPM",
  "Viral_percent_non_mosquito",
  "Viral_RPM_non_mosquito",
  "Average_coverage",
  "Normalized_mosquito_coverage",
  "log10_Viral_Reads",
  "log10_Other_Reads",
  "log10_Non_mosquito_Reads",
  "log10_Mosquito_Reads",
  "log10_Total_Reads",
  "log10_Viral_RPM",
  "log10_Other_RPM",
  "log10_Viral_RPM_non_mosquito"
)

missing_metrics <- setdiff(metrics_to_analyze, names(df_paired))

if (length(missing_metrics) > 0) {
  stop(
    "The following metrics are missing from df_paired: ",
    paste(missing_metrics, collapse = ", ")
  )
}

long_df <- df_paired %>%
  dplyr::select(Pair_ID, Sample, Host, Enrichment, dplyr::all_of(metrics_to_analyze)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(metrics_to_analyze),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  dplyr::arrange(Metric, Pair_ID, Enrichment)

# ------------------------------------------------------------
# Calculate summary statistics
# ------------------------------------------------------------

summary_stats <- long_df %>%
  dplyr::group_by(Metric, Enrichment) %>%
  dplyr::summarise(
    n = sum(!is.na(Value)),
    mean = mean(Value, na.rm = TRUE),
    median = median(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    min = min(Value, na.rm = TRUE),
    max = max(Value, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  summary_stats,
  file.path(tables_dir, "summary_statistics_by_enrichment.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Run paired Wilcoxon tests and effect-size calculations
# ------------------------------------------------------------

paired_tests <- long_df %>%
  dplyr::filter(!is.na(Value)) %>%
  dplyr::arrange(Metric, Pair_ID, Enrichment) %>%
  dplyr::group_by(Metric) %>%
  rstatix::wilcox_test(Value ~ Enrichment, paired = TRUE) %>%
  dplyr::ungroup() %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  tibble::as_tibble()

effect_sizes <- long_df %>%
  dplyr::filter(!is.na(Value)) %>%
  dplyr::arrange(Metric, Pair_ID, Enrichment) %>%
  dplyr::group_by(Metric) %>%
  rstatix::wilcox_effsize(Value ~ Enrichment, paired = TRUE) %>%
  dplyr::ungroup() %>%
  tibble::as_tibble()

paired_results <- paired_tests %>%
  dplyr::left_join(
    effect_sizes %>%
      dplyr::select(Metric, effsize, magnitude),
    by = "Metric"
  ) %>%
  dplyr::rename(
    wilcoxon_p = p,
    wilcoxon_p_adj_BH = p.adj,
    wilcoxon_p_adj_signif = p.adj.signif,
    effect_size_r = effsize,
    effect_size_magnitude = magnitude
  )

write.csv(
  paired_results,
  file.path(tables_dir, "paired_wilcoxon_results.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Calculate paired change caused by enrichment
# ------------------------------------------------------------

paired_change_wide <- long_df %>%
  dplyr::select(Pair_ID, Host, Enrichment, Metric, Value) %>%
  dplyr::mutate(
    Enrichment = as.character(Enrichment)
  ) %>%
  tidyr::pivot_wider(
    id_cols = c(Pair_ID, Host, Metric),
    names_from = Enrichment,
    values_from = Value
  )

if (!"FALSE" %in% names(paired_change_wide)) {
  paired_change_wide[["FALSE"]] <- NA_real_
}

if (!"TRUE" %in% names(paired_change_wide)) {
  paired_change_wide[["TRUE"]] <- NA_real_
}

paired_change_wide <- paired_change_wide %>%
  dplyr::rename(
    unenriched_value = `FALSE`,
    enriched_value = `TRUE`
  ) %>%
  dplyr::mutate(
    absolute_change = enriched_value - unenriched_value,
    percent_change = safe_percent_change(enriched_value, unenriched_value),
    fold_change = safe_fold_change(enriched_value, unenriched_value),
    log2_fold_change = dplyr::case_when(
      is.na(fold_change) ~ NA_real_,
      fold_change <= 0 ~ NA_real_,
      TRUE ~ log2(fold_change)
    )
  )

write.csv(
  paired_change_wide,
  file.path(tables_dir, "paired_enrichment_change_by_sample.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Average enrichment-associated change across all paired samples
# ------------------------------------------------------------

paired_change_summary_overall <- paired_change_wide %>%
  dplyr::group_by(Metric) %>%
  dplyr::summarise(
    n_pairs = dplyr::n(),
    n_pairs_with_percent_change = sum(!is.na(percent_change)),
    mean_unenriched = mean(unenriched_value, na.rm = TRUE),
    mean_enriched = mean(enriched_value, na.rm = TRUE),
    median_unenriched = median(unenriched_value, na.rm = TRUE),
    median_enriched = median(enriched_value, na.rm = TRUE),
    mean_absolute_change = mean(absolute_change, na.rm = TRUE),
    median_absolute_change = median(absolute_change, na.rm = TRUE),
    mean_percent_change = mean(percent_change, na.rm = TRUE),
    median_percent_change = median(percent_change, na.rm = TRUE),
    mean_fold_change = mean(fold_change, na.rm = TRUE),
    median_fold_change = median(fold_change, na.rm = TRUE),
    mean_log2_fold_change = mean(log2_fold_change, na.rm = TRUE),
    median_log2_fold_change = median(log2_fold_change, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    paired_results %>%
      dplyr::select(
        Metric,
        wilcoxon_p,
        wilcoxon_p_adj_BH,
        wilcoxon_p_adj_signif,
        effect_size_r,
        effect_size_magnitude
      ),
    by = "Metric"
  ) %>%
  dplyr::arrange(Metric)

write.csv(
  paired_change_summary_overall,
  file.path(tables_dir, "paired_enrichment_change_summary_overall.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Average enrichment-associated change by mosquito host
# ------------------------------------------------------------

paired_change_summary_by_host <- paired_change_wide %>%
  dplyr::group_by(Host, Metric) %>%
  dplyr::summarise(
    n_pairs = dplyr::n(),
    n_pairs_with_percent_change = sum(!is.na(percent_change)),
    mean_unenriched = mean(unenriched_value, na.rm = TRUE),
    mean_enriched = mean(enriched_value, na.rm = TRUE),
    median_unenriched = median(unenriched_value, na.rm = TRUE),
    median_enriched = median(enriched_value, na.rm = TRUE),
    mean_absolute_change = mean(absolute_change, na.rm = TRUE),
    median_absolute_change = median(absolute_change, na.rm = TRUE),
    mean_percent_change = mean(percent_change, na.rm = TRUE),
    median_percent_change = median(percent_change, na.rm = TRUE),
    mean_fold_change = mean(fold_change, na.rm = TRUE),
    median_fold_change = median(fold_change, na.rm = TRUE),
    mean_log2_fold_change = mean(log2_fold_change, na.rm = TRUE),
    median_log2_fold_change = median(log2_fold_change, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Metric, Host)

write.csv(
  paired_change_summary_by_host,
  file.path(tables_dir, "paired_enrichment_change_summary_by_host.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Main read-change summary for manuscript-relevant metrics
# ------------------------------------------------------------

main_change_metrics <- c(
  "Viral_reads",
  "Viral_RPM",
  "Percent_viral",
  "Viral_RPM_non_mosquito",
  "Non_mosquito_reads",
  "Non_mosquito_RPM",
  "Other_reads",
  "Other_RPM",
  "Mosquito_reads",
  "Mosquito_RPM",
  "Percent_mosquito",
  "Average_coverage",
  "Normalized_mosquito_coverage",
  "Total_reads"
)

main_change_summary <- paired_change_summary_overall %>%
  dplyr::filter(Metric %in% main_change_metrics) %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      mean_absolute_change > 0 ~ "Higher after enrichment",
      mean_absolute_change < 0 ~ "Lower after enrichment",
      mean_absolute_change == 0 ~ "No average change",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::arrange(match(Metric, main_change_metrics))

write.csv(
  main_change_summary,
  file.path(tables_dir, "main_enrichment_average_change_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Excel workbook for change summaries
# ------------------------------------------------------------

change_sheets <- list(
  change_by_sample = paired_change_wide,
  overall_summary = paired_change_summary_overall,
  summary_by_host = paired_change_summary_by_host,
  main_metrics = main_change_summary
)

wb_change <- openxlsx::createWorkbook()

for (sheet_name in names(change_sheets)) {
  openxlsx::addWorksheet(wb_change, sheet_name)
  openxlsx::writeData(wb_change, sheet_name, change_sheets[[sheet_name]])
}

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom"
)

for (sheet_name in names(change_sheets)) {
  n_cols <- ncol(change_sheets[[sheet_name]])
  if (n_cols > 0) {
    openxlsx::addStyle(
      wb_change,
      sheet = sheet_name,
      style = header_style,
      rows = 1,
      cols = 1:n_cols,
      gridExpand = TRUE
    )
    openxlsx::freezePane(wb_change, sheet = sheet_name, firstRow = TRUE)
    openxlsx::setColWidths(wb_change, sheet = sheet_name, cols = 1:n_cols, widths = "auto")
  }
}

openxlsx::saveWorkbook(
  wb_change,
  file.path(tables_dir, "paired_enrichment_average_change_summary.xlsx"),
  overwrite = TRUE
)

# ------------------------------------------------------------
# Build main interpretation table
# ------------------------------------------------------------

main_metrics <- c(
  "Viral_reads",
  "Percent_viral",
  "Viral_RPM",
  "log10_Viral_Reads",
  "Viral_percent_non_mosquito",
  "Viral_RPM_non_mosquito",
  "Average_coverage",
  "Normalized_mosquito_coverage",
  "Total_reads",
  "Non_mosquito_reads",
  "Other_reads",
  "Other_RPM",
  "Percent_other",
  "Percent_mosquito",
  "Mosquito_reads"
)

main_interpretation_table <- summary_stats %>%
  dplyr::filter(Metric %in% main_metrics) %>%
  dplyr::select(Metric, Enrichment, median, mean, sd) %>%
  tidyr::pivot_wider(
    names_from = Enrichment,
    values_from = c(median, mean, sd),
    names_sep = "_"
  ) %>%
  dplyr::left_join(
    paired_results %>%
      dplyr::select(
        Metric,
        wilcoxon_p,
        wilcoxon_p_adj_BH,
        wilcoxon_p_adj_signif,
        effect_size_r,
        effect_size_magnitude
      ),
    by = "Metric"
  ) %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      wilcoxon_p_adj_signif != "ns" & median_TRUE > median_FALSE ~ "Higher in enriched",
      wilcoxon_p_adj_signif != "ns" & median_TRUE < median_FALSE ~ "Lower in enriched",
      TRUE ~ "Not significant"
    ),
    interpretation_note = dplyr::case_when(
      Metric %in% c("Average_coverage", "Normalized_mosquito_coverage") ~
        "Mosquito genome coverage metric; do not interpret as viral genome coverage.",
      Metric %in% c(
        "Viral_reads",
        "Percent_viral",
        "Viral_RPM",
        "log10_Viral_Reads",
        "Viral_percent_non_mosquito",
        "Viral_RPM_non_mosquito"
      ) ~
        "Viral read recovery metric.",
      Metric %in% c("Other_reads", "Other_RPM", "Percent_other") ~
        "Non-mosquito, non-viral read metric.",
      Metric %in% c("Percent_mosquito", "Mosquito_reads") ~
        "Mosquito/host read metric.",
      TRUE ~
        "Read-composition metric."
    )
  ) %>%
  dplyr::arrange(match(Metric, main_metrics))

write.csv(
  main_interpretation_table,
  file.path(tables_dir, "main_interpretation_table.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Define plotting colors and themes
# ------------------------------------------------------------

enrichment_colors <- c(
  "FALSE" = "red",
  "TRUE" = "steelblue"
)

host_colors <- c(
  "Aedes albopictus" = "#E41A1C",
  "Culex pipiens" = "#377EB8",
  "Culiseta melanura" = "#4DAF4A"
)

base_enrichment_theme <- ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 20, hjust = 0, margin = ggplot2::margin(b = 8)),
    axis.title = ggplot2::element_text(size = 17),
    axis.text = ggplot2::element_text(size = 14),
    strip.background = ggplot2::element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.8
    ),
    strip.text = ggplot2::element_text(size = 13, face = "bold"),
    legend.position = "none",
    panel.spacing = grid::unit(1.4, "lines"),
    plot.margin = ggplot2::margin(t = 12, r = 24, b = 12, l = 12)
  )

base_host_legend_theme <- base_enrichment_theme +
  ggplot2::theme(
    legend.position = "bottom",
    legend.title = ggplot2::element_text(size = 13, face = "bold"),
    legend.text = ggplot2::element_text(size = 12)
  )

# ------------------------------------------------------------
# Plot read composition and mosquito genome coverage
# ------------------------------------------------------------

composition_metrics <- c(
  "Normalized_mosquito_coverage",
  "Other_RPM",
  "Percent_mosquito",
  "Viral_RPM"
)

composition_metric_labels <- c(
  Normalized_mosquito_coverage = "Mosquito Genome Coverage\nper Million Mosquito Reads",
  Other_RPM = "Other Reads\nper Million Total Reads",
  Percent_mosquito = "Percent Mosquito Reads",
  Viral_RPM = "Viral Reads\nper Million Total Reads"
)

composition_df <- long_df %>%
  dplyr::filter(Metric %in% composition_metrics) %>%
  dplyr::mutate(
    Metric_label = dplyr::recode(Metric, !!!composition_metric_labels)
  )

composition_stats <- paired_results %>%
  dplyr::filter(Metric %in% composition_metrics) %>%
  dplyr::mutate(
    Metric_label = dplyr::recode(Metric, !!!composition_metric_labels),
    stat_label = format_stat_label(wilcoxon_p_adj_signif, effect_size_r)
  )

composition_stats <- add_stat_positions(composition_df, composition_stats)

p_composition_box <- ggplot2::ggplot(
  composition_df,
  ggplot2::aes(x = Enrichment, y = Value, color = Enrichment)
) +
  ggplot2::geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    linewidth = 0.8
  ) +
  ggplot2::geom_jitter(
    width = 0.10,
    size = 2.6,
    alpha = 0.75
  ) +
  ggplot2::facet_wrap(
    ~ Metric_label,
    scales = "free_y",
    ncol = 2
  ) +
  ggplot2::geom_text(
    data = composition_stats,
    ggplot2::aes(
      x = 1.5,
      y = y_pos,
      label = stat_label
    ),
    inherit.aes = FALSE,
    size = 4.8,
    fontface = "bold"
  ) +
  ggplot2::scale_color_manual(values = enrichment_colors) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.28))) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    title = "Effect of Enrichment on Read Composition\nand Mosquito Genome Coverage",
    x = "Enrichment",
    y = "Measurement"
  ) +
  base_enrichment_theme

ggplot2::ggsave(
  filename = file.path(figures_dir, "enrichment_boxplot_read_composition_host_coverage.png"),
  plot = p_composition_box,
  width = 11,
  height = 8,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(figures_dir, "enrichment_boxplot_read_composition_host_coverage.pdf"),
  plot = p_composition_box,
  width = 11,
  height = 8
)

p_composition_paired <- ggplot2::ggplot(
  composition_df,
  ggplot2::aes(
    x = Enrichment,
    y = Value,
    group = Pair_ID
  )
) +
  ggplot2::geom_line(
    color = "gray65",
    linewidth = 0.6,
    alpha = 0.65
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = Host),
    size = 3.2,
    alpha = 0.9
  ) +
  ggplot2::facet_wrap(
    ~ Metric_label,
    scales = "free_y",
    ncol = 2
  ) +
  ggplot2::geom_text(
    data = composition_stats,
    ggplot2::aes(
      x = 1.5,
      y = y_pos,
      label = stat_label
    ),
    inherit.aes = FALSE,
    size = 4.8,
    fontface = "bold",
    color = "black"
  ) +
  ggplot2::scale_color_manual(values = host_colors, drop = FALSE) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.28))) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    title = "Effect of Enrichment on Read Composition\nand Mosquito Genome Coverage",
    x = "Enrichment",
    y = "Measurement",
    color = "Mosquito host"
  ) +
  base_host_legend_theme

ggplot2::ggsave(
  filename = file.path(figures_dir, "enrichment_paired_read_composition_host_coverage_by_host.png"),
  plot = p_composition_paired,
  width = 11,
  height = 8,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(figures_dir, "enrichment_paired_read_composition_host_coverage_by_host.pdf"),
  plot = p_composition_paired,
  width = 11,
  height = 8
)

ggplot2::ggsave(
  filename = file.path(figures_dir, "enrichment_paired_read_composition_host_coverage.png"),
  plot = p_composition_paired,
  width = 11,
  height = 8,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(figures_dir, "enrichment_paired_read_composition_host_coverage.pdf"),
  plot = p_composition_paired,
  width = 11,
  height = 8
)

# ------------------------------------------------------------
# Plot viral read recovery
# ------------------------------------------------------------

viral_metrics <- c(
  "Percent_viral",
  "Viral_RPM",
  "log10_Viral_Reads"
)

viral_metric_labels <- c(
  Percent_viral = "Percent Viral Reads",
  Viral_RPM = "Viral Reads per Million",
  log10_Viral_Reads = "Viral Reads (log10)"
)

viral_plot_df <- long_df %>%
  dplyr::filter(Metric %in% viral_metrics) %>%
  dplyr::mutate(
    Metric_label = dplyr::recode(Metric, !!!viral_metric_labels)
  )

viral_plot_stats <- paired_results %>%
  dplyr::filter(Metric %in% viral_metrics) %>%
  dplyr::mutate(
    Metric_label = dplyr::recode(Metric, !!!viral_metric_labels),
    stat_label = format_stat_label(wilcoxon_p_adj_signif, effect_size_r)
  )

viral_plot_stats <- add_stat_positions(viral_plot_df, viral_plot_stats)

p_viral <- ggplot2::ggplot(
  viral_plot_df,
  ggplot2::aes(
    x = Enrichment,
    y = Value,
    group = Pair_ID
  )
) +
  ggplot2::geom_line(
    color = "gray65",
    linewidth = 0.6,
    alpha = 0.65
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = Host),
    size = 3.2,
    alpha = 0.9
  ) +
  ggplot2::facet_wrap(
    ~ Metric_label,
    scales = "free_y",
    ncol = 2
  ) +
  ggplot2::geom_text(
    data = viral_plot_stats,
    ggplot2::aes(
      x = 1.5,
      y = y_pos,
      label = stat_label
    ),
    inherit.aes = FALSE,
    size = 4.8,
    fontface = "bold",
    color = "black"
  ) +
  ggplot2::scale_color_manual(values = host_colors, drop = FALSE) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.28))) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    title = "Effect of Enrichment on Viral Read Recovery",
    x = "Enrichment",
    y = "Measurement",
    color = "Mosquito host"
  ) +
  base_host_legend_theme

ggplot2::ggsave(
  filename = file.path(figures_dir, "enrichment_paired_viral_read_recovery.png"),
  plot = p_viral,
  width = 11,
  height = 6.5,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(figures_dir, "enrichment_paired_viral_read_recovery.pdf"),
  plot = p_viral,
  width = 11,
  height = 6.5
)

# ------------------------------------------------------------
# Plot mosquito genome coverage
# ------------------------------------------------------------

host_coverage_metrics <- c(
  "Average_coverage",
  "Normalized_mosquito_coverage"
)

host_coverage_metric_labels <- c(
  Average_coverage = "Mosquito Genome Coverage\n(Average Depth)",
  Normalized_mosquito_coverage = "Mosquito Genome Coverage\nper Million Mosquito Reads"
)

host_coverage_df <- long_df %>%
  dplyr::filter(Metric %in% host_coverage_metrics) %>%
  dplyr::mutate(
    Metric_label = dplyr::recode(Metric, !!!host_coverage_metric_labels)
  )

host_coverage_stats <- paired_results %>%
  dplyr::filter(Metric %in% host_coverage_metrics) %>%
  dplyr::mutate(
    Metric_label = dplyr::recode(Metric, !!!host_coverage_metric_labels),
    stat_label = format_stat_label(wilcoxon_p_adj_signif, effect_size_r)
  )

host_coverage_stats <- add_stat_positions(host_coverage_df, host_coverage_stats)

p_host_coverage <- ggplot2::ggplot(
  host_coverage_df,
  ggplot2::aes(
    x = Enrichment,
    y = Value,
    group = Pair_ID
  )
) +
  ggplot2::geom_line(
    color = "gray65",
    linewidth = 0.6,
    alpha = 0.65
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = Host),
    size = 3.2,
    alpha = 0.9
  ) +
  ggplot2::facet_wrap(
    ~ Metric_label,
    scales = "free_y",
    ncol = 2
  ) +
  ggplot2::geom_text(
    data = host_coverage_stats,
    ggplot2::aes(
      x = 1.5,
      y = y_pos,
      label = stat_label
    ),
    inherit.aes = FALSE,
    size = 4.8,
    fontface = "bold",
    color = "black"
  ) +
  ggplot2::scale_color_manual(values = host_colors, drop = FALSE) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.28))) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    title = "Effect of Enrichment on Mosquito Genome Coverage",
    x = "Enrichment",
    y = "Measurement",
    color = "Mosquito host"
  ) +
  base_host_legend_theme

ggplot2::ggsave(
  filename = file.path(figures_dir, "enrichment_paired_mosquito_genome_coverage.png"),
  plot = p_host_coverage,
  width = 11,
  height = 5.5,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(figures_dir, "enrichment_paired_mosquito_genome_coverage.pdf"),
  plot = p_host_coverage,
  width = 11,
  height = 5.5
)

# ------------------------------------------------------------
# Plot paired enrichment effects by mosquito host
# ------------------------------------------------------------

host_species_metrics <- c(
  "Percent_viral",
  "Viral_RPM",
  "Percent_mosquito",
  "Normalized_mosquito_coverage",
  "Other_RPM"
)

host_species_labels <- c(
  Percent_viral = "Percent Viral Reads",
  Viral_RPM = "Viral Reads\nper Million Total Reads",
  Percent_mosquito = "Percent Mosquito Reads",
  Normalized_mosquito_coverage = "Mosquito Genome Coverage\nper Million Mosquito Reads",
  Other_RPM = "Other Reads\nper Million Total Reads"
)

host_species_df <- long_df %>%
  dplyr::filter(Metric %in% host_species_metrics) %>%
  dplyr::mutate(
    Metric_label = dplyr::recode(Metric, !!!host_species_labels),
    Host = factor(
      Host,
      levels = c(
        "Aedes albopictus",
        "Culex pipiens",
        "Culiseta melanura"
      )
    ),
    Host_Metric = paste(Host, Metric_label, sep = "\n")
  )

p_host_species <- ggplot2::ggplot(
  host_species_df,
  ggplot2::aes(
    x = Enrichment,
    y = Value,
    group = Pair_ID
  )
) +
  ggplot2::geom_line(
    color = "gray65",
    linewidth = 0.6,
    alpha = 0.65
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = Host),
    size = 3.0,
    alpha = 0.9
  ) +
  ggplot2::facet_wrap(
    ~ Host_Metric,
    scales = "free_y",
    ncol = 5
  ) +
  ggplot2::scale_color_manual(values = host_colors, drop = FALSE) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.20))) +
  ggplot2::labs(
    title = "Paired Enrichment Effects by Mosquito Host",
    x = "Enrichment",
    y = "Measurement",
    color = "Mosquito host"
  ) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 20, hjust = 0, margin = ggplot2::margin(b = 8)),
    axis.title = ggplot2::element_text(size = 17),
    axis.text = ggplot2::element_text(size = 13),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    strip.background = ggplot2::element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.8
    ),
    strip.text = ggplot2::element_text(size = 8.5, face = "bold"),
    panel.spacing = grid::unit(1.0, "lines"),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(size = 13),
    legend.text = ggplot2::element_text(size = 12),
    plot.margin = ggplot2::margin(t = 12, r = 18, b = 12, l = 12)
  )

ggplot2::ggsave(
  filename = file.path(figures_dir, "enrichment_paired_by_host_species.png"),
  plot = p_host_species,
  width = 18,
  height = 10,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(figures_dir, "enrichment_paired_by_host_species.pdf"),
  plot = p_host_species,
  width = 18,
  height = 10
)

# ------------------------------------------------------------
# Save run information
# ------------------------------------------------------------

sink(file.path(logs_dir, "session_info.txt"))
cat("Enrichment analysis completed:\n")
cat(as.character(Sys.time()), "\n\n")

cat("Input file:\n")
cat(normalizePath(input_file), "\n\n")

cat("Output folder:\n")
cat(normalizePath(output_dir), "\n\n")

cat("Complete pairs used:\n")
print(complete_pairs)

cat("\n\nColumn names used:\n")
print(names(df_clean))

cat("\n\nImportant interpretation note:\n")
cat("Viral_RPM = viral reads per million total reads.\n")
cat("Other_RPM = other reads per million total reads.\n")
cat("Average_coverage and Normalized_mosquito_coverage are mosquito genome coverage metrics.\n")
cat("Normalized_mosquito_coverage = average mosquito genome coverage per million mosquito-mapped reads.\n")
cat("Coverage metrics should not be interpreted as viral genome coverage or viral recovery metrics.\n")
cat("Percent/fold change values are NA when the unenriched value was zero and enriched value was greater than zero.\n")

cat("\n\nSession info:\n")
sessionInfo()
sink()

# ------------------------------------------------------------
# Print completion message
# ------------------------------------------------------------

cat("\nAnalysis complete.\n")
cat("Results saved to:\n")
cat(normalizePath(output_dir), "\n\n")

cat("Main output tables:\n")
cat(file.path(tables_dir, "paired_wilcoxon_results.csv"), "\n")
cat(file.path(tables_dir, "main_interpretation_table.csv"), "\n")
cat(file.path(tables_dir, "summary_statistics_by_enrichment.csv"), "\n")
cat(file.path(tables_dir, "pair_check.csv"), "\n")
cat(file.path(tables_dir, "paired_enrichment_change_by_sample.csv"), "\n")
cat(file.path(tables_dir, "paired_enrichment_change_summary_overall.csv"), "\n")
cat(file.path(tables_dir, "paired_enrichment_change_summary_by_host.csv"), "\n")
cat(file.path(tables_dir, "main_enrichment_average_change_summary.csv"), "\n")
cat(file.path(tables_dir, "paired_enrichment_average_change_summary.xlsx"), "\n\n")

cat("Main output figures:\n")
cat(file.path(figures_dir, "enrichment_boxplot_read_composition_host_coverage.png"), "\n")
cat(file.path(figures_dir, "enrichment_paired_read_composition_host_coverage_by_host.png"), "\n")
cat(file.path(figures_dir, "enrichment_paired_viral_read_recovery.png"), "\n")
cat(file.path(figures_dir, "enrichment_paired_mosquito_genome_coverage.png"), "\n")
cat(file.path(figures_dir, "enrichment_paired_by_host_species.png"), "\n")