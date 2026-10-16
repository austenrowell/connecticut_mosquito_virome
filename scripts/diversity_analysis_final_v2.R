# ----------------------------
# OPTIONAL ENTREZ API KEY SETUP
# ----------------------------
# Recommended: set your NCBI API key for efficiency
# ----------------------------
# USER SETTINGS
# ----------------------------

alb_dir  <- "ALB_annotated/"
cpip_dir <- "CPIP_annotated/"
cmel_dir <- "CMEL_annotated/"

metadata_file <- "sample_to_host_clean.csv"
collection_metadata_file <- "sample_collection_metadata.csv"
summary_file  <- "summary_complete.csv"
outdir <- "diversity_results/"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# OUTPUT DIRECTORY STRUCTURE
# ----------------------------

make_analysis_dirs <- function(analysis_name) {
  base_dir <- file.path(outdir, analysis_name)
  dir_list <- list(
    base = base_dir,
    data = file.path(base_dir, "data"),
    results = file.path(base_dir, "results"),
    figures = file.path(base_dir, "figures")
  )

  purrr::walk(dir_list, ~dir.create(.x, showWarnings = FALSE, recursive = TRUE))
  dir_list
}

abundance_dirs <- make_analysis_dirs("01_abundance_matrices")
alpha_dirs     <- make_analysis_dirs("02_alpha_diversity")
pcoa_dirs      <- make_analysis_dirs("03_pcoa_ordination")
beta_dirs      <- make_analysis_dirs("04_beta_diversity_tests")
permdisp_dirs  <- make_analysis_dirs("05_permdisp")
heatmap_dirs   <- make_analysis_dirs("06_heatmaps")
deseq_dirs     <- make_analysis_dirs("07_differential_abundance")
family_dirs    <- make_analysis_dirs("08_family_endogenization")
coverage_dirs  <- make_analysis_dirs("09_coverage_analysis")
spacetime_dirs <- make_analysis_dirs("10_spatiotemporal_analysis")
logs_dirs      <- make_analysis_dirs("99_logs")

# ----------------------------
# PACKAGES
# ----------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ape)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(DESeq2)
  library(taxize)
  library(forcats)
  library(scales)
  library(grid)
  library(broom)
  library(ggpubr)
  library(rstatix)
  library(lubridate)
})

# ----------------------------
# CHECK ENTREZ API KEY
# ----------------------------

entrez_key <- Sys.getenv("ENTREZ_KEY", unset = "")

if (nzchar(entrez_key)) {
  message("ENTREZ_KEY detected in environment. Taxonomy lookup will use the configured key.")
} else {
  message("No ENTREZ_KEY detected. Taxonomy lookup will run without a key and may be slower/rate-limited.")
}

# ----------------------------
# LOAD METADATA
# ----------------------------

meta <- readr::read_csv(metadata_file, show_col_types = FALSE) |>
  dplyr::rename(sample = Sample, host = Host)

if (file.exists(collection_metadata_file)) {
  collection_meta <- readr::read_csv(collection_metadata_file, show_col_types = FALSE) |>
    dplyr::rename(
      sample = Sample,
      collection_host = Host,
      town = Town,
      site = Site,
      latitude = Latitude,
      longitude = Longitude,
      collection_date_raw = Date,
      county = County,
      trap_type = Trap_type,
      mosquitoes = Mosquitoes
    ) |>
    dplyr::mutate(
      collection_date = suppressWarnings(as.Date(lubridate::parse_date_time(
        collection_date_raw,
        orders = c("ymd", "mdy", "dmy", "dbY", "bdY", "d-b-Y", "m/d/Y")
      ))),
      collection_year = lubridate::year(collection_date),
      collection_month_num = lubridate::month(collection_date),
      collection_month = factor(
        lubridate::month(collection_date, label = TRUE, abbr = FALSE),
        levels = month.name,
        ordered = TRUE
      ),
      latitude = as.numeric(latitude),
      longitude = as.numeric(longitude),
      site = as.factor(site),
      county = as.factor(county)
    ) |>
    dplyr::select(
      sample,
      town,
      site,
      county,
      latitude,
      longitude,
      collection_date,
      collection_year,
      collection_month_num,
      collection_month,
      trap_type,
      mosquitoes
    )

  meta <- meta |>
    dplyr::left_join(collection_meta, by = "sample")

  readr::write_tsv(
    meta,
    file.path(spacetime_dirs$data, "metadata_with_collection_site_county_month.tsv")
  )
} else {
  warning("collection_metadata_file was not found. Spatial-temporal analyses will be skipped.")
}

host_dir_map <- tibble::tribble(
  ~host, ~dir,
  "Aedes albopictus",  alb_dir,
  "Culex pipiens",     cpip_dir,
  "Culiseta melanura", cmel_dir
)

meta <- meta |>
  dplyr::left_join(host_dir_map, by = "host")

if (any(is.na(meta$dir))) {
  stop("Some hosts in metadata do not match known host_dir_map.")
}

# ----------------------------
# FUNCTION: Parse Kraken report
# ----------------------------

parse_kreport <- function(file) {
  sample <- gsub("\\.k2report$", "", basename(file))
  df <- readr::read_tsv(file, col_names = FALSE, show_col_types = FALSE)
  
  colnames(df) <- c(
    "pct", "clade_reads", "taxon_reads",
    "clade_frags", "taxon_frags",
    "rank", "taxid", "name"
  )
  
  df |>
    dplyr::filter(stringr::str_detect(rank, "^S")) |>
    dplyr::mutate(
      sample = sample,
      endo = stringr::str_detect(name, "\\*$"),
      taxon = stringr::str_trim(stringr::str_remove(name, "\\*$"))
    ) |>
    dplyr::select(sample, taxon, clade_reads, endo)
}

parse_kreport_with_taxid <- function(file) {
  sample <- gsub("\\.k2report$", "", basename(file))
  df <- readr::read_tsv(file, col_names = FALSE, show_col_types = FALSE)
  
  colnames(df) <- c(
    "pct", "clade_reads", "taxon_reads",
    "clade_frags", "taxon_frags",
    "rank", "taxid", "name"
  )
  
  df |>
    dplyr::filter(stringr::str_detect(rank, "^S")) |>
    dplyr::mutate(
      sample = sample,
      taxid = as.character(taxid),
      endo = stringr::str_detect(name, "\\*$"),
      taxon = stringr::str_trim(stringr::str_remove(name, "\\*$"))
    ) |>
    dplyr::select(sample, taxid, taxon, clade_reads, endo)
}

# ----------------------------
# FIND VALID FILES
# ----------------------------

valid_files <- meta |>
  dplyr::mutate(path = file.path(dir, paste0(sample, ".k2report"))) |>
  dplyr::filter(file.exists(path)) |>
  dplyr::group_by(sample) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup()

cat("Loading", nrow(valid_files), "validated files\n")

if (nrow(valid_files) == 0) {
  stop("No valid Kraken2 report files were found.")
}

# ----------------------------
# LOAD ALL DATA
# ----------------------------

long <- purrr::map_dfr(valid_files$path, parse_kreport)
readr::write_tsv(long, file.path(abundance_dirs$data, "abundance_long.tsv"))

# ----------------------------
# BUILD MATRICES
# ----------------------------

make_matrix <- function(df) {
  df |>
    dplyr::group_by(sample, taxon) |>
    dplyr::summarise(clade_reads = sum(clade_reads), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = taxon,
      values_from = clade_reads,
      values_fill = 0
    ) |>
    tibble::column_to_rownames("sample")
}

total_mat   <- make_matrix(long)
endo_mat    <- make_matrix(dplyr::filter(long, endo))
nonendo_mat <- make_matrix(dplyr::filter(long, !endo))

readr::write_tsv(as.data.frame(total_mat),   file.path(abundance_dirs$data, "abundance_total.tsv"))
readr::write_tsv(as.data.frame(endo_mat),    file.path(abundance_dirs$data, "abundance_endo.tsv"))
readr::write_tsv(as.data.frame(nonendo_mat), file.path(abundance_dirs$data, "abundance_nonendo.tsv"))

# ----------------------------
# ALPHA DIVERSITY
# ----------------------------

alpha_div <- function(mat) {
  tibble::tibble(
    sample   = rownames(mat),
    richness = vegan::specnumber(mat),
    shannon  = vegan::diversity(mat, index = "shannon"),
    simpson  = vegan::diversity(mat, index = "simpson")
  )
}

alpha_total   <- alpha_div(total_mat)
alpha_endo    <- alpha_div(endo_mat)
alpha_nonendo <- alpha_div(nonendo_mat)

readr::write_tsv(alpha_total,   file.path(alpha_dirs$data, "alpha_total.tsv"))
readr::write_tsv(alpha_endo,    file.path(alpha_dirs$data, "alpha_endo.tsv"))
readr::write_tsv(alpha_nonendo, file.path(alpha_dirs$data, "alpha_nonendo.tsv"))

# ----------------------------
# PREP META FOR PLOTTING
# ----------------------------

meta_for_plot <- meta |>
  dplyr::distinct(
    sample,
    host,
    town,
    site,
    county,
    latitude,
    longitude,
    collection_date,
    collection_year,
    collection_month_num,
    collection_month,
    trap_type,
    mosquitoes
  )

add_host <- function(mat, meta_tbl) {
  df <- as.data.frame(mat)
  df$sample <- rownames(df)
  dplyr::left_join(df, meta_tbl, by = "sample")
}

total_df   <- add_host(total_mat, meta_for_plot)
endo_df    <- add_host(endo_mat, meta_for_plot)
nonendo_df <- add_host(nonendo_mat, meta_for_plot)

# ----------------------------
# PCoA + VARIANCE + ENVFIT
# ----------------------------

run_pcoa_color <- function(df, label) {
  metadata_cols <- c(
    "sample",
    "host",
    "town",
    "site",
    "county",
    "trap_type",
    "mosquitoes",
    "date",
    "Date",
    "collection_date",
    "collection_date_raw",
    "collection_year",
    "collection_month",
    "collection_month_num",
    "latitude",
    "longitude",
    "Latitude",
    "Longitude",
    "dir",
    "path"
  )

  mat_df <- df |>
    dplyr::select(-dplyr::any_of(metadata_cols))

  non_numeric_cols <- names(mat_df)[!vapply(mat_df, is.numeric, logical(1))]
  if (length(non_numeric_cols) > 0) {
    message(
      "Removing non-numeric columns from the ", label,
      " abundance matrix before Bray-Curtis analysis: ",
      paste(non_numeric_cols, collapse = ", ")
    )
    mat_df <- mat_df |>
      dplyr::select(where(is.numeric))
  }

  if (ncol(mat_df) == 0) {
    stop("No numeric abundance columns remain for ", label, " after removing metadata columns.")
  }

  mat <- as.matrix(mat_df)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- df$sample

  mat[is.na(mat)] <- 0

  zero_rows <- which(rowSums(mat, na.rm = TRUE) == 0)
  if (length(zero_rows) > 0) {
    mat[zero_rows, ] <- mat[zero_rows, ] + 1e-6
  }

  dist <- vegan::vegdist(mat, method = "bray")
  attr(dist, "Labels") <- rownames(mat)

  pcoa_res <- ape::pcoa(dist)

  available_axes <- min(5, ncol(pcoa_res$vectors))
  ord <- pcoa_res$vectors[, seq_len(available_axes), drop = FALSE]
  var_explained <- pcoa_res$values$Relative_eig[seq_len(available_axes)]

  readr::write_tsv(
    tibble::tibble(
      axis = paste0("PCoA", seq_len(available_axes)),
      variance_explained = var_explained
    ),
    file.path(pcoa_dirs$data, paste0("pcoa_", label, "_variance.tsv"))
  )

  plot_df <- as.data.frame(ord)
  colnames(plot_df) <- paste0("PCoA", seq_len(ncol(plot_df)))
  plot_df$sample <- rownames(mat)
  plot_df$host <- df$host

  readr::write_tsv(
    plot_df,
    file.path(pcoa_dirs$data, paste0("pcoa_", label, "_host.tsv"))
  )

  meta_model <- df |>
    dplyr::select(sample, host)
  rownames(meta_model) <- meta_model$sample

  ef_axes <- as.matrix(ord[, 1:min(2, ncol(ord)), drop = FALSE])
  ef <- vegan::envfit(ef_axes ~ host, data = meta_model, permutations = 999)
  ef_df <- as.data.frame(ef$factors$pvals)
  colnames(ef_df) <- "p_value"
  ef_df$variable <- rownames(ef_df)

  readr::write_tsv(
    ef_df,
    file.path(pcoa_dirs$results, paste0("pcoa_", label, "_envfit.tsv"))
  )

  host_colors <- c(
    "Aedes albopictus" = "#E41A1C",
    "Culex pipiens" = "#377EB8",
    "Culiseta melanura" = "#4DAF4A"
  )

  make_pcoa_plot <- function(plot_data, x_axis, y_axis, labeled = FALSE) {
    x_num <- as.integer(stringr::str_remove(x_axis, "PCoA"))
    y_num <- as.integer(stringr::str_remove(y_axis, "PCoA"))

    p <- ggplot(plot_data, aes(x = .data[[x_axis]], y = .data[[y_axis]], color = host)) +
      geom_point(size = ifelse(labeled, 3, 3.5), alpha = 0.9) +
      theme_minimal(base_size = 14) +
      labs(
        title = paste("Bray-Curtis PCoA:", label),
        x = paste0(x_axis, " (", round(var_explained[x_num] * 100, 1), "%)"),
        y = paste0(y_axis, " (", round(var_explained[y_num] * 100, 1), "%)"),
        color = "Host"
      ) +
      scale_color_manual(values = host_colors, drop = FALSE) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "right",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )

    if (labeled) {
      p <- p + ggrepel::geom_text_repel(aes(label = sample), size = 3, max.overlaps = Inf)
    }

    p
  }

  axis_pairs <- list(
    c("PCoA1", "PCoA2"),
    c("PCoA1", "PCoA3"),
    c("PCoA2", "PCoA3")
  )

  for (axis_pair in axis_pairs) {
    x_axis <- axis_pair[1]
    y_axis <- axis_pair[2]

    if (!all(c(x_axis, y_axis) %in% colnames(plot_df))) {
      next
    }

    axis_label <- paste0(x_axis, "_", y_axis)

    p_labeled <- make_pcoa_plot(plot_df, x_axis, y_axis, labeled = TRUE)
    p_unlabeled <- make_pcoa_plot(plot_df, x_axis, y_axis, labeled = FALSE)

    if (axis_label == "PCoA1_PCoA2") {
      ggsave(
        file.path(pcoa_dirs$figures, paste0("pcoa_", label, "_host.pdf")),
        p_labeled,
        width = 7,
        height = 5,
        bg = "white"
      )
    }

    ggsave(
      file.path(pcoa_dirs$figures, paste0("pcoa_", label, "_host_", axis_label, "_labeled.pdf")),
      p_labeled,
      width = 7,
      height = 5,
      bg = "white"
    )

    ggsave(
      file.path(pcoa_dirs$figures, paste0("pcoa_", label, "_host_", axis_label, "_labeled.png")),
      p_labeled,
      width = 7,
      height = 5,
      dpi = 300,
      bg = "white"
    )

    ggsave(
      file.path(pcoa_dirs$figures, paste0("pcoa_", label, "_host_", axis_label, "_unlabeled.pdf")),
      p_unlabeled,
      width = 7,
      height = 5,
      bg = "white"
    )

    ggsave(
      file.path(pcoa_dirs$figures, paste0("pcoa_", label, "_host_", axis_label, "_unlabeled.png")),
      p_unlabeled,
      width = 7,
      height = 5,
      dpi = 300,
      bg = "white"
    )

    if (requireNamespace("svglite", quietly = TRUE)) {
      ggsave(
        file.path(pcoa_dirs$figures, paste0("pcoa_", label, "_host_", axis_label, "_unlabeled.svg")),
        p_unlabeled,
        width = 7,
        height = 5,
        device = svglite::svglite,
        bg = "white"
      )
    }
  }

  if (all(c("PCoA1", "PCoA2", "PCoA3") %in% colnames(plot_df)) &&
      requireNamespace("plotly", quietly = TRUE) &&
      requireNamespace("htmlwidgets", quietly = TRUE)) {

    p3d <- plotly::plot_ly(
      data = plot_df,
      x = ~PCoA1,
      y = ~PCoA2,
      z = ~PCoA3,
      color = ~host,
      colors = host_colors,
      text = ~sample,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 5)
    ) |>
      plotly::layout(
        title = paste("3D Bray-Curtis PCoA:", label),
        scene = list(
          xaxis = list(title = paste0("PCoA1 (", round(var_explained[1] * 100, 1), "%)")),
          yaxis = list(title = paste0("PCoA2 (", round(var_explained[2] * 100, 1), "%)")),
          zaxis = list(title = paste0("PCoA3 (", round(var_explained[3] * 100, 1), "%)"))
        ),
        legend = list(title = list(text = "Host"))
      )

    htmlwidgets::saveWidget(
      p3d,
      file.path(pcoa_dirs$figures, paste0("pcoa_", label, "_3D_interactive.html")),
      selfcontained = TRUE
    )
  }

  dist
}

dist_total   <- run_pcoa_color(total_df, "total")
dist_endo    <- run_pcoa_color(endo_df, "endo")
dist_nonendo <- run_pcoa_color(nonendo_df, "nonendo")


# ----------------------------
# PERMANOVA + ANOSIM
# ----------------------------

run_tests <- function(dist, meta_tbl, label) {
  df <- meta_tbl[match(labels(dist), meta_tbl$sample), , drop = FALSE]
  df$host <- factor(df$host)
  
  perm <- vegan::adonis2(dist ~ host, data = df)
  capture.output(perm, file = file.path(beta_dirs$results, paste0("permanova_", label, ".txt")))
  
  ano <- vegan::anosim(dist, df$host)
  capture.output(ano, file = file.path(beta_dirs$results, paste0("anosim_", label, ".txt")))
}

run_tests(dist_total, meta_for_plot, "total")
run_tests(dist_endo, meta_for_plot, "endo")
run_tests(dist_nonendo, meta_for_plot, "nonendo")

# ----------------------------
# PAIRWISE PERMANOVA BY HOST
# ----------------------------

run_pairwise_permanova <- function(dist_obj, meta_tbl, label, permutations = 999) {
  
  df <- meta_tbl |>
    dplyr::filter(sample %in% labels(dist_obj)) |>
    dplyr::distinct(sample, host)
  
  df <- df[match(labels(dist_obj), df$sample), , drop = FALSE]
  df$host <- factor(df$host)
  
  host_levels <- levels(df$host)
  host_pairs <- combn(host_levels, 2, simplify = FALSE)
  
  pairwise_results <- purrr::map_dfr(host_pairs, function(pair) {
    
    host_1 <- pair[1]
    host_2 <- pair[2]
    
    keep_samples <- df |>
      dplyr::filter(host %in% c(host_1, host_2)) |>
      dplyr::pull(sample)
    
    pair_df <- df |>
      dplyr::filter(sample %in% keep_samples) |>
      dplyr::mutate(host = droplevels(host))
    
    pair_dist <- as.dist(as.matrix(dist_obj)[keep_samples, keep_samples])
    
    perm <- vegan::adonis2(
      pair_dist ~ host,
      data = pair_df,
      permutations = permutations
    )
    
    perm_df <- as.data.frame(perm)
    
    r2_col <- if ("R2" %in% colnames(perm_df)) "R2" else NA_character_
    f_col <- if ("F" %in% colnames(perm_df)) "F" else NA_character_
    p_col <- grep("Pr", colnames(perm_df), value = TRUE)[1]
    
    result_row <- tibble::tibble(
      dataset = label,
      comparison = paste(host_1, "vs", host_2),
      host_1 = host_1,
      host_2 = host_2,
      n_host_1 = sum(pair_df$host == host_1),
      n_host_2 = sum(pair_df$host == host_2),
      r_squared = if (!is.na(r2_col)) perm_df[1, r2_col] else NA_real_,
      f_statistic = if (!is.na(f_col)) perm_df[1, f_col] else NA_real_,
      p_value = if (!is.na(p_col)) perm_df[1, p_col] else NA_real_
    )
    
    safe_label <- paste0(
      "pairwise_permanova_",
      label,
      "_",
      gsub(" ", "_", host_1),
      "_vs_",
      gsub(" ", "_", host_2),
      ".txt"
    )
    
    capture.output(
      perm,
      file = file.path(beta_dirs$results, safe_label)
    )
    
    result_row
  })
  
  pairwise_results <- pairwise_results |>
    dplyr::mutate(
      p_adj_BH = p.adjust(p_value, method = "BH"),
      significance = dplyr::case_when(
        is.na(p_adj_BH) ~ NA_character_,
        p_adj_BH <= 0.001 ~ "***",
        p_adj_BH <= 0.01 ~ "**",
        p_adj_BH <= 0.05 ~ "*",
        p_adj_BH <= 0.1 ~ ".",
        TRUE ~ "ns"
      )
    )
  
  readr::write_tsv(
    pairwise_results,
    file.path(beta_dirs$results, paste0("pairwise_permanova_", label, "_summary.tsv"))
  )
  
  pairwise_results
}

pairwise_permanova_total <- run_pairwise_permanova(
  dist_total,
  meta_for_plot,
  "total"
)

pairwise_permanova_endo <- run_pairwise_permanova(
  dist_endo,
  meta_for_plot,
  "endo"
)

pairwise_permanova_nonendo <- run_pairwise_permanova(
  dist_nonendo,
  meta_for_plot,
  "nonendo"
)

pairwise_permanova_all <- dplyr::bind_rows(
  pairwise_permanova_total,
  pairwise_permanova_endo,
  pairwise_permanova_nonendo
)

readr::write_tsv(
  pairwise_permanova_all,
  file.path(beta_dirs$results, "pairwise_permanova_all_summary.tsv")
)

cat("\nPairwise PERMANOVA complete.\n")
cat("Results written to:", beta_dirs$results, "\n")

# ----------------------------
# SPATIAL-TEMPORAL PERMANOVA + ENVFIT
# Date/location variables are used only here and in PCoA-based outputs.
# ----------------------------

subset_dist_by_samples <- function(dist_obj, samples) {
  dist_mat <- as.matrix(dist_obj)
  keep <- intersect(labels(dist_obj), samples)
  as.dist(dist_mat[keep, keep, drop = FALSE])
}

run_spatiotemporal_tests <- function(dist_obj, meta_tbl, label) {
  required_cols <- c(
    "sample",
    "host",
    "site",
    "county",
    "latitude",
    "longitude",
    "collection_month_num",
    "collection_month"
  )

  missing_cols <- setdiff(required_cols, colnames(meta_tbl))
  if (length(missing_cols) > 0) {
    message("Spatial-temporal tests skipped for ", label, ": missing columns: ", paste(missing_cols, collapse = ", "))
    return(NULL)
  }

  df <- meta_tbl |>
    dplyr::filter(sample %in% labels(dist_obj)) |>
    dplyr::mutate(
      host = factor(
        host,
        levels = c(
          "Aedes albopictus",
          "Culex pipiens",
          "Culiseta melanura"
        )
      ),
      site = factor(site),
      county = factor(county),
      collection_month = factor(collection_month, levels = month.name, ordered = TRUE),
      collection_month_num = as.numeric(collection_month_num),
      latitude = as.numeric(latitude),
      longitude = as.numeric(longitude)
    ) |>
    dplyr::filter(
      !is.na(host),
      !is.na(latitude),
      !is.na(longitude),
      !is.na(collection_month_num),
      !is.na(collection_month),
      !is.na(county),
      !is.na(site)
    )

  if (nrow(df) < 6) {
    message("Spatial-temporal tests skipped for ", label, ": too few complete samples.")
    return(NULL)
  }

  df <- df[match(intersect(labels(dist_obj), df$sample), df$sample), , drop = FALSE]
  dist_sub <- subset_dist_by_samples(dist_obj, df$sample)

  continuous_model <- vegan::adonis2(
    dist_sub ~ host + latitude + longitude + collection_month_num,
    data = df,
    by = "margin",
    permutations = 999
  )

  county_month_model <- vegan::adonis2(
    dist_sub ~ host + county + collection_month,
    data = df,
    by = "margin",
    permutations = 999
  )

  site_month_model <- vegan::adonis2(
    dist_sub ~ host + site + collection_month,
    data = df,
    by = "margin",
    permutations = 999
  )

  capture.output(
    continuous_model,
    file = file.path(spacetime_dirs$results, paste0("permanova_", label, "_host_lat_long_month.txt"))
  )

  capture.output(
    county_month_model,
    file = file.path(spacetime_dirs$results, paste0("permanova_", label, "_host_county_month.txt"))
  )

  capture.output(
    site_month_model,
    file = file.path(spacetime_dirs$results, paste0("permanova_", label, "_host_site_month.txt"))
  )

  model_summary <- dplyr::bind_rows(
    broom::tidy(continuous_model) |>
      dplyr::mutate(model = "host + latitude + longitude + month_numeric"),
    broom::tidy(county_month_model) |>
      dplyr::mutate(model = "host + county + month"),
    broom::tidy(site_month_model) |>
      dplyr::mutate(model = "host + site + month")
  ) |>
    dplyr::mutate(dataset = label) |>
    dplyr::select(dataset, model, dplyr::everything())

  readr::write_tsv(
    model_summary,
    file.path(spacetime_dirs$results, paste0("permanova_", label, "_spatiotemporal_summary.tsv"))
  )

  pcoa_res <- ape::pcoa(dist_sub)
  ord <- as.data.frame(pcoa_res$vectors[, 1:min(3, ncol(pcoa_res$vectors)), drop = FALSE])
  colnames(ord) <- paste0("PCoA", seq_len(ncol(ord)))
  ord$sample <- rownames(ord)

  plot_df <- ord |>
    dplyr::left_join(df, by = "sample") |>
    dplyr::mutate(
      county = factor(county),
      site = factor(site),
      collection_month = factor(collection_month, levels = month.name, ordered = TRUE),
      host = factor(
        host,
        levels = c(
          "Aedes albopictus",
          "Culex pipiens",
          "Culiseta melanura"
        )
      )
    )

  readr::write_tsv(
    plot_df,
    file.path(spacetime_dirs$data, paste0("pcoa_", label, "_spatiotemporal_scores.tsv"))
  )

  var_explained <- pcoa_res$values$Relative_eig[1:min(3, ncol(pcoa_res$vectors))]

  ef_numeric <- vegan::envfit(
    pcoa_res$vectors[, 1:2, drop = FALSE] ~ latitude + longitude + collection_month_num,
    data = df,
    permutations = 999
  )

  ef_factors <- vegan::envfit(
    pcoa_res$vectors[, 1:2, drop = FALSE] ~ county + site + collection_month,
    data = df,
    permutations = 999
  )

  capture.output(
    ef_numeric,
    file = file.path(spacetime_dirs$results, paste0("envfit_", label, "_numeric_location_month.txt"))
  )

  capture.output(
    ef_factors,
    file = file.path(spacetime_dirs$results, paste0("envfit_", label, "_factor_site_county_month.txt"))
  )

  host_shapes <- c(
    "Aedes albopictus" = 16,
    "Culex pipiens" = 17,
    "Culiseta melanura" = 15
  )

  make_spatiotemporal_pcoa_plot <- function(plot_data, x_axis, y_axis, color_var, labeled = FALSE) {
    x_num <- as.integer(stringr::str_remove(x_axis, "PCoA"))
    y_num <- as.integer(stringr::str_remove(y_axis, "PCoA"))

    color_label <- dplyr::case_when(
      color_var == "county" ~ "County",
      color_var == "site" ~ "Site",
      color_var == "collection_month" ~ "Collection month",
      TRUE ~ color_var
    )

    p <- ggplot(
      plot_data,
      aes(
        x = .data[[x_axis]],
        y = .data[[y_axis]],
        color = .data[[color_var]],
        shape = host
      )
    ) +
      geom_point(size = 3.2, alpha = 0.9) +
      scale_shape_manual(values = host_shapes, drop = FALSE) +
      theme_minimal(base_size = 14) +
      labs(
        title = paste("Bray-Curtis PCoA:", label, "colored by", color_label),
        x = paste0(x_axis, " (", round(var_explained[x_num] * 100, 1), "%)"),
        y = paste0(y_axis, " (", round(var_explained[y_num] * 100, 1), "%)"),
        color = color_label,
        shape = "Host"
      ) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "right",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )

    if (labeled) {
      p <- p +
        ggrepel::geom_text_repel(
          aes(label = sample),
          size = 3,
          max.overlaps = Inf,
          show.legend = FALSE
        )
    }

    p
  }

  axis_pairs <- list(
    c("PCoA1", "PCoA2"),
    c("PCoA1", "PCoA3"),
    c("PCoA2", "PCoA3")
  )

  spatiotemporal_vars <- c("county", "site", "collection_month")

  for (axis_pair in axis_pairs) {
    x_axis <- axis_pair[1]
    y_axis <- axis_pair[2]

    if (!all(c(x_axis, y_axis) %in% colnames(plot_df))) {
      next
    }

    axis_label <- paste0(x_axis, "_", y_axis)

    for (sp_var in spatiotemporal_vars) {
      if (!sp_var %in% colnames(plot_df)) {
        next
      }

      if (all(is.na(plot_df[[sp_var]]))) {
        next
      }

      safe_var <- stringr::str_replace_all(sp_var, "[^A-Za-z0-9_]+", "_")

      p_labeled <- make_spatiotemporal_pcoa_plot(
        plot_data = plot_df,
        x_axis = x_axis,
        y_axis = y_axis,
        color_var = sp_var,
        labeled = TRUE
      )

      p_unlabeled <- make_spatiotemporal_pcoa_plot(
        plot_data = plot_df,
        x_axis = x_axis,
        y_axis = y_axis,
        color_var = sp_var,
        labeled = FALSE
      )

      ggsave(
        file.path(
          spacetime_dirs$figures,
          paste0("pcoa_", label, "_colored_by_", safe_var, "_shape_by_host_", axis_label, "_labeled.pdf")
        ),
        p_labeled,
        width = 7,
        height = 5,
        bg = "white"
      )

      ggsave(
        file.path(
          spacetime_dirs$figures,
          paste0("pcoa_", label, "_colored_by_", safe_var, "_shape_by_host_", axis_label, "_labeled.png")
        ),
        p_labeled,
        width = 7,
        height = 5,
        dpi = 300,
        bg = "white"
      )

      ggsave(
        file.path(
          spacetime_dirs$figures,
          paste0("pcoa_", label, "_colored_by_", safe_var, "_shape_by_host_", axis_label, "_unlabeled.pdf")
        ),
        p_unlabeled,
        width = 7,
        height = 5,
        bg = "white"
      )

      ggsave(
        file.path(
          spacetime_dirs$figures,
          paste0("pcoa_", label, "_colored_by_", safe_var, "_shape_by_host_", axis_label, "_unlabeled.png")
        ),
        p_unlabeled,
        width = 7,
        height = 5,
        dpi = 300,
        bg = "white"
      )

      if (requireNamespace("svglite", quietly = TRUE)) {
        ggsave(
          file.path(
            spacetime_dirs$figures,
            paste0("pcoa_", label, "_colored_by_", safe_var, "_shape_by_host_", axis_label, "_unlabeled.svg")
          ),
          p_unlabeled,
          width = 7,
          height = 5,
          device = svglite::svglite,
          bg = "white"
        )
      }
    }
  }

  # Interactive 3D spatiotemporal plots are intentionally skipped.
  # The 2D spatiotemporal PCoA outputs retain host identifiers using point shape,
  # while color represents county, site, or collection month.

  readr::write_tsv(
    df,
    file.path(spacetime_dirs$data, paste0("metadata_used_for_spatiotemporal_", label, ".tsv"))
  )

  invisible(
    list(
      continuous_model = continuous_model,
      county_month_model = county_month_model,
      site_month_model = site_month_model,
      envfit_numeric = ef_numeric,
      envfit_factors = ef_factors,
      metadata = df,
      pcoa_scores = plot_df
    )
  )
}

spacetime_total <- run_spatiotemporal_tests(dist_total, meta_for_plot, "total")
spacetime_endo <- run_spatiotemporal_tests(dist_endo, meta_for_plot, "endo")
spacetime_nonendo <- run_spatiotemporal_tests(dist_nonendo, meta_for_plot, "nonendo")

# ----------------------------
# ENHANCED PERMDISP
# ----------------------------

calc_effect_sizes_from_betadisper_anova <- function(anova_table) {
  ss_between <- anova_table[1, "Sum Sq"]
  ss_resid   <- anova_table[2, "Sum Sq"]
  df_between <- anova_table[1, "Df"]
  df_resid   <- anova_table[2, "Df"]
  
  ss_total <- ss_between + ss_resid
  eta2 <- ss_between / ss_total
  ms_res <- ss_resid / df_resid
  omega2 <- (ss_between - df_between * ms_res) / (ss_total + ms_res)
  if (omega2 < 0) omega2 <- 0
  
  list(
    eta2 = as.numeric(eta2),
    omega2 = as.numeric(omega2)
  )
}

run_permdisp_enhanced <- function(dist, meta_tbl, label) {
  df <- meta_tbl[match(labels(dist), meta_tbl$sample), , drop = FALSE]
  df$host <- factor(df$host)
  
  if (length(unique(df$host)) < 2) {
    message("PERMDISP skipped for ", label, ": fewer than 2 groups")
    return(NULL)
  }
  
  bd <- vegan::betadisper(dist, df$host)
  perm <- vegan::permutest(bd, permutations = 999)
  
  anova_bd <- anova(bd)
  anova_df <- as.data.frame(anova_bd)
  
  if (!"Sum Sq" %in% colnames(anova_df)) {
    idx <- colnames(anova_df) %in% c("SumSq", "Sumsq")
    colnames(anova_df)[idx] <- "Sum Sq"
  }
  
  eff <- calc_effect_sizes_from_betadisper_anova(anova_df)
  
  dispersion_df <- data.frame(
    sample = names(bd$distances),
    distance = bd$distances,
    host = df$host
  )
  
  group_summary <- dispersion_df |>
    dplyr::group_by(host) |>
    dplyr::summarise(
      mean_distance = mean(distance),
      sd_distance = sd(distance),
      n = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(desc(mean_distance))
  
  tuk <- tryCatch(TukeyHSD(bd), error = function(e) NULL)
  pair_t <- tryCatch(
    pairwise.t.test(dispersion_df$distance, dispersion_df$host, p.adjust.method = "BH", pool.sd = FALSE),
    error = function(e) NULL
  )
  pair_w <- tryCatch(
    pairwise.wilcox.test(dispersion_df$distance, dispersion_df$host, p.adjust.method = "BH"),
    error = function(e) NULL
  )
  
  results_prefix <- file.path(permdisp_dirs$results, paste0("permdisp_", label))
  figure_prefix <- file.path(permdisp_dirs$figures, paste0("permdisp_", label))
  data_prefix <- file.path(permdisp_dirs$data, paste0("permdisp_", label))
  capture.output(perm, file = paste0(results_prefix, ".txt"))
  readr::write_tsv(group_summary, paste0(data_prefix, "_group_means.tsv"))
  capture.output(anova_bd, file = paste0(results_prefix, "_anova.txt"))
  if (!is.null(tuk)) capture.output(tuk, file = paste0(results_prefix, "_tukey.txt"))
  if (!is.null(pair_t)) capture.output(pair_t, file = paste0(results_prefix, "_pairwise_t_test.txt"))
  if (!is.null(pair_w)) capture.output(pair_w, file = paste0(results_prefix, "_pairwise_wilcox.txt"))
  
  p <- ggplot(dispersion_df, aes(x = host, y = distance, color = host)) +
    geom_jitter(width = 0.15, size = 2.5) +
    geom_boxplot(alpha = 0.25, outlier.shape = NA) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Distance to Centroid:", label),
      x = "Host",
      y = "Distance to centroid"
    ) +
    scale_color_brewer(palette = "Set1")
  
  ggsave(
    paste0(figure_prefix, "_distance_boxplot.pdf"),
    p,
    width = 7,
    height = 5,
    bg = "white"
  )
  
  pdf(paste0(figure_prefix, "_betadisper_plot.pdf"), width = 7, height = 5)
  plot(bd, main = paste("PERMDISP:", label))
  dev.off()
  
  top_host_row <- group_summary[1, ]
  top_host <- as.character(top_host_row$host)
  top_mean <- top_host_row$mean_distance
  
  significant_pairs <- character(0)
  
  if (!is.null(pair_t) && !is.null(pair_t$p.value)) {
    pmat <- pair_t$p.value
    which_sig <- which(pmat < 0.05, arr.ind = TRUE)
    if (nrow(which_sig) > 0) {
      for (i in seq_len(nrow(which_sig))) {
        r <- which_sig[i, 1]
        c <- which_sig[i, 2]
        host1 <- rownames(pmat)[r]
        host2 <- colnames(pmat)[c]
        pair <- paste(host1, "vs", host2)
        if (grepl(top_host, pair, fixed = TRUE)) {
          significant_pairs <- c(significant_pairs, pair)
        }
      }
    }
  }
  
  if (length(significant_pairs) == 0 && !is.null(tuk)) {
    tuk_tab <- tuk[[1]]
    if (!is.null(tuk_tab) && "p adj" %in% colnames(tuk_tab)) {
      sig_idx <- which(tuk_tab[, "p adj"] < 0.05)
      if (length(sig_idx) > 0) {
        comps <- rownames(tuk_tab)[sig_idx]
        significant_pairs <- comps[grepl(top_host, comps, fixed = TRUE)]
      }
    }
  }
  
  result <- list(
    label = label,
    permutest = perm,
    anova = anova_bd,
    effect_sizes = eff,
    group_summary = group_summary,
    tukey = tuk,
    pairwise_t_test = pair_t,
    pairwise_wilcox = pair_w,
    dispersion_df = dispersion_df,
    top_host = top_host,
    top_mean = top_mean,
    significant_pairs_for_top = unique(significant_pairs)
  )
  
  summary_lines <- c(
    paste0("Dataset: ", label),
    paste0("PERMDISP permutations summary written to: ", paste0(results_prefix, ".txt")),
    paste0("ANOVA table written to: ", paste0(results_prefix, "_anova.txt")),
    paste0("Effect sizes: eta2 = ", round(eff$eta2, 4), ", omega2 = ", round(eff$omega2, 4)),
    paste0("Top host by mean distance: ", top_host, " (mean distance = ", round(top_mean, 4), ")"),
    "Significant pairwise comparisons involving top host (BH-corrected t-test or Tukey):",
    if (length(result$significant_pairs_for_top) == 0) {
      "  - none (no significant pairwise comparisons involving top host at p < 0.05)"
    } else {
      paste0("  - ", result$significant_pairs_for_top)
    }
  )
  
  writeLines(summary_lines, con = paste0(results_prefix, "_summary.txt"))
  
  result
}

res_total   <- run_permdisp_enhanced(dist_total, meta_for_plot, "total")
res_endo    <- run_permdisp_enhanced(dist_endo, meta_for_plot, "endo")
res_nonendo <- run_permdisp_enhanced(dist_nonendo, meta_for_plot, "nonendo")

# ----------------------------
# COLLATE PERMDISP SUMMARY
# ----------------------------

collect_summary_row <- function(res) {
  if (is.null(res)) return(NULL)
  
  tibble::tibble(
    dataset = res$label,
    perm_p_value = {
      pval <- NA_real_
      try({
        if (!is.null(res$permutest) && !is.null(res$permutest$tab)) {
          tab <- as.data.frame(res$permutest$tab)
          if ("Pr(>F)" %in% colnames(tab)) {
            pval <- as.numeric(tab[1, "Pr(>F)"])
          } else if ("Pr..F." %in% colnames(tab)) {
            pval <- as.numeric(tab[1, "Pr..F."])
          } else if (ncol(tab) >= 4) {
            pval <- as.numeric(tab[1, ncol(tab)])
          }
        }
      }, silent = TRUE)
      pval
    },
    anova_p_value = {
      pval_a <- NA_real_
      try({
        adf <- as.data.frame(res$anova)
        pcol <- grep("Pr", colnames(adf), value = TRUE)
        if (length(pcol) >= 1) pval_a <- as.numeric(adf[1, pcol[1]])
      }, silent = TRUE)
      pval_a
    },
    eta2 = res$effect_sizes$eta2,
    omega2 = res$effect_sizes$omega2,
    top_host = res$top_host,
    top_host_mean_distance = res$top_mean,
    n_groups = nrow(res$group_summary),
    significant_pairs_involving_top = paste(res$significant_pairs_for_top, collapse = "; ")
  )
}

summary_rows <- dplyr::bind_rows(
  collect_summary_row(res_total),
  collect_summary_row(res_endo),
  collect_summary_row(res_nonendo)
)

readr::write_tsv(summary_rows, file.path(permdisp_dirs$results, "permdisp_overview_summary.tsv"))

# ----------------------------
# HEATMAPS
# ----------------------------

prop_mat <- total_mat
row_sums <- rowSums(prop_mat)
prop_mat[row_sums == 0, ] <- prop_mat[row_sums == 0, ] + 1e-6
prop_mat <- prop_mat / rowSums(prop_mat)

bc_dist <- as.matrix(vegan::vegdist(prop_mat, method = "bray"))
readr::write_tsv(as.data.frame(bc_dist), file.path(heatmap_dirs$data, "bray_curtis_distance_matrix.tsv"))

# Heatmap graphics use host only.
# Date, site, county, latitude, and longitude are reserved for PCoA/spatiotemporal statistics.
heatmap_annotation <- meta |>
  dplyr::filter(sample %in% rownames(bc_dist)) |>
  dplyr::select(sample, Host = host) |>
  dplyr::mutate(Host = factor(Host))

annotation_df <- heatmap_annotation |>
  tibble::column_to_rownames("sample")

sample_order <- heatmap_annotation |>
  dplyr::arrange(Host, sample) |>
  dplyr::pull(sample)

bc_ordered <- bc_dist[sample_order, sample_order, drop = FALSE]
ann_ordered <- annotation_df[sample_order, , drop = FALSE]

heat_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")))(100)

pheatmap::pheatmap(
  bc_ordered,
  filename = file.path(heatmap_dirs$figures, "heatmap_bray_ordered_by_host.pdf"),
  color = heat_colors,
  annotation_row = ann_ordered,
  annotation_col = ann_ordered,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6,
  fontsize_col = 6,
  angle_col = 45,
  width = 12,
  height = 12,
  main = "Bray-Curtis distance (ordered by host)"
)

pheatmap::pheatmap(
  bc_dist,
  filename = file.path(heatmap_dirs$figures, "heatmap_bray_clustered_by_host_annotation.pdf"),
  color = heat_colors,
  annotation_row = annotation_df,
  annotation_col = annotation_df,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6,
  fontsize_col = 6,
  angle_col = 45,
  width = 12,
  height = 12,
  main = "Bray-Curtis distance (clustered)"
)

# ----------------------------
# BOXPLOTS: ALPHA DIVERSITY BY HOST
# ----------------------------

alpha_total_host <- dplyr::left_join(alpha_total, meta_for_plot, by = "sample")
alpha_endo_host  <- dplyr::left_join(alpha_endo,  meta_for_plot, by = "sample")

plot_alpha <- function(df, metric, label) {
  p <- ggplot(df, aes(x = host, y = .data[[metric]], color = host)) +
    geom_jitter(width = 0.15, size = 2.5) +
    geom_boxplot(alpha = 0.25, outlier.shape = NA) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste(label, metric, "by host"),
      x = "Host",
      y = metric
    ) +
    scale_color_brewer(palette = "Set1")
  
  ggsave(
    file.path(alpha_dirs$figures, paste0("boxplot_", label, "_", metric, "_by_host.pdf")),
    plot = p,
    width = 7,
    height = 5,
    bg = "white"
  )
}

plot_alpha(alpha_total_host, "shannon",  "total")
plot_alpha(alpha_total_host, "richness", "total")
plot_alpha(alpha_endo_host,  "shannon",  "endo")
plot_alpha(alpha_endo_host,  "richness", "endo")

# ----------------------------
# DIFFERENTIAL ABUNDANCE (DESeq2)
# ----------------------------

counts <- as.matrix(total_mat)
counts <- t(counts)
storage.mode(counts) <- "integer"

keep <- rowSums(counts >= 10) >= 3
counts <- counts[keep, , drop = FALSE]

meta_deseq <- meta_for_plot |>
  dplyr::filter(sample %in% colnames(counts)) |>
  tibble::column_to_rownames("sample")

meta_deseq <- meta_deseq[colnames(counts), , drop = FALSE]

host_levels_original <- unique(meta_deseq$host)
host_levels_safe <- make.names(host_levels_original)
host_map <- stats::setNames(host_levels_original, host_levels_safe)

meta_deseq$host_safe <- factor(make.names(meta_deseq$host), levels = host_levels_safe)

stopifnot(identical(colnames(counts), rownames(meta_deseq)))

dds_lrt <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta_deseq,
  design = ~ host_safe
)

dds_lrt <- DESeq2::DESeq(dds_lrt, test = "LRT", reduced = ~ 1, fitType = "local")

res_lrt <- DESeq2::results(dds_lrt)
res_lrt_df <- as.data.frame(res_lrt) |>
  tibble::rownames_to_column("taxon")

readr::write_tsv(res_lrt_df, file.path(deseq_dirs$results, "deseq_host_effect_all_taxa.tsv"))

dds_wald <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta_deseq,
  design = ~ host_safe
)

dds_wald <- DESeq2::DESeq(dds_wald, fitType = "local")

make_volcano <- function(res_df, label,
                         padj_cutoff = 0.05,
                         lfc_cutoff = 1,
                         top_n = 20,
                         outdir = outdir) {
  
  if (!"taxon" %in% colnames(res_df)) {
    if (!is.null(rownames(res_df))) {
      res_df <- tibble::rownames_to_column(res_df, "taxon")
    } else {
      stop("res_df must have rownames or a 'taxon' column")
    }
  }
  
  if (!all(c("log2FoldChange", "padj") %in% colnames(res_df))) {
    stop("res_df must contain 'log2FoldChange' and 'padj' columns")
  }
  
  nonzero_padj <- res_df$padj[!is.na(res_df$padj) & res_df$padj > 0]
  min_nonzero <- if (length(nonzero_padj) > 0) min(nonzero_padj) else 1e-300
  floor_val <- max(min_nonzero * 0.1, 1e-300)
  
  res_df <- res_df |>
    dplyr::mutate(
      padj_plot = dplyr::case_when(
        is.na(padj) ~ 1,
        padj == 0 ~ floor_val,
        TRUE ~ padj
      ),
      negLog10Padj = -log10(padj_plot),
      absLFC = abs(log2FoldChange),
      is_sig = ifelse(!is.na(padj) & padj < padj_cutoff & absLFC >= lfc_cutoff, "yes", "no")
    )
  
  finite_values <- res_df$negLog10Padj[is.finite(res_df$negLog10Padj)]
  if (length(finite_values) > 0) {
    cap_val <- max(10, stats::quantile(finite_values, 0.99, na.rm = TRUE) * 2)
    res_df$negLog10Padj <- pmin(res_df$negLog10Padj, cap_val)
  }
  
  top_taxa <- res_df |>
    dplyr::filter(is_sig == "yes") |>
    dplyr::arrange(desc(absLFC)) |>
    dplyr::slice_head(n = top_n)
  
  xlim_max <- max(
    3,
    stats::quantile(abs(res_df$log2FoldChange[is.finite(res_df$log2FoldChange)]), 0.99, na.rm = TRUE) * 1.1
  )
  xlim <- c(-xlim_max, xlim_max)
  
  ylim_max <- max(1, max(res_df$negLog10Padj[is.finite(res_df$negLog10Padj)], na.rm = TRUE) * 1.05)
  ylim <- c(0, ylim_max)
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = negLog10Padj)) +
    geom_point(aes(color = is_sig), alpha = 0.8, size = 1.6, show.legend = FALSE) +
    scale_color_manual(values = c("no" = "grey70", "yes" = "red")) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "blue") +
    {
      if (nrow(top_taxa) > 0) {
        ggrepel::geom_text_repel(
          data = top_taxa,
          aes(label = taxon),
          size = 3,
          max.overlaps = Inf,
          box.padding = 0.3,
          segment.alpha = 0.6
        )
      } else {
        NULL
      }
    } +
    theme_minimal(base_size = 13) +
    labs(
      title = paste("Volcano:", label),
      x = "log2 fold change",
      y = "-log10(adjusted p-value)"
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))
  
  out_pdf <- file.path(deseq_dirs$figures, paste0("volcano_", label, ".pdf"))
  out_png <- file.path(deseq_dirs$figures, paste0("volcano_", label, ".png"))
  
  ggsave(out_pdf, p, width = 8, height = 6, bg = "white")
  ggsave(out_png, p, width = 8, height = 6, dpi = 200, bg = "white")
  
  invisible(list(plot = p, top_taxa = top_taxa, pdf = out_pdf, png = out_png))
}

make_all_volcanoes <- function(dds, outdir, host_map, padj_cutoff = 0.05, lfc_cutoff = 1, top_n = 20) {
  hosts_safe <- levels(SummarizedExperiment::colData(dds)$host_safe)
  
  if (length(hosts_safe) < 2) {
    stop("Need at least two host groups to generate contrasts.")
  }
  
  message("Generating volcano plots for all host contrasts:")
  
  combs <- combn(hosts_safe, 2, simplify = FALSE)
  
  for (pair in combs) {
    h1_safe <- pair[1]
    h2_safe <- pair[2]
    h1 <- unname(host_map[h1_safe])
    h2 <- unname(host_map[h2_safe])
    
    label <- paste0(gsub(" ", "_", h1), "_vs_", gsub(" ", "_", h2))
    message("  - ", label)
    
    res <- DESeq2::results(dds, contrast = c("host_safe", h1_safe, h2_safe))
    
    res_df <- as.data.frame(res) |>
      tibble::rownames_to_column("taxon")
    
    readr::write_tsv(res_df, file.path(deseq_dirs$results, paste0("deseq_", label, ".tsv")))
    
    make_volcano(
      res_df,
      label,
      padj_cutoff = padj_cutoff,
      lfc_cutoff = lfc_cutoff,
      top_n = top_n,
      outdir = outdir
    )
  }
  
  message("All volcano plots and DESeq2 tables saved to: ", outdir)
}

res_culex_vs_aedes <- DESeq2::results(
  dds_wald,
  contrast = c("host_safe", make.names("Culex pipiens"), make.names("Aedes albopictus"))
)

res_culiseta_vs_aedes <- DESeq2::results(
  dds_wald,
  contrast = c("host_safe", make.names("Culiseta melanura"), make.names("Aedes albopictus"))
)

sig <- subset(res_culex_vs_aedes, padj < 0.05)
sig$enriched_in <- ifelse(sig$log2FoldChange > 0, "Culex pipiens", "Aedes albopictus")

sig1 <- subset(res_culiseta_vs_aedes, padj < 0.05)
sig1$enriched_in <- ifelse(sig1$log2FoldChange > 0, "Culiseta melanura", "Aedes albopictus")

print(res_culex_vs_aedes)
print(res_culiseta_vs_aedes)

make_all_volcanoes(
  dds_wald,
  outdir,
  host_map = host_map,
  padj_cutoff = 0.05,
  lfc_cutoff = 1,
  top_n = 20
)

# ----------------------------
# FAMILY-LEVEL ENDOGENIZATION COMPARISON
# ----------------------------

cat("Re-parsing Kraken reports to obtain taxids...\n")

long_tax <- purrr::map_dfr(valid_files$path, parse_kreport_with_taxid)

readr::write_tsv(long_tax, file.path(family_dirs$data, "abundance_long_with_taxid.tsv"))

entrez_key <- Sys.getenv("ENTREZ_KEY", unset = "")

if (nzchar(entrez_key)) {
  message("Using ENTREZ_KEY from environment.")
} else {
  message("No ENTREZ_KEY found in environment; taxonomy lookup may be slower and more rate-limited.")
  message("Set it in ~/.Renviron as: ENTREZ_KEY=your_ncbi_api_key, then restart R.")
}

taxid_vec <- unique(long_tax$taxid)
taxid_vec <- as.character(taxid_vec[!is.na(taxid_vec) & taxid_vec != ""])

cat("Resolving families for", length(taxid_vec), "taxids...\n")

get_family_safe <- function(taxids, chunk_size = 50, sleep_sec = 0.4, retry_singletons = TRUE) {
  taxids <- unique(as.character(taxids))
  taxids <- taxids[!is.na(taxids) & taxids != ""]
  
  fam <- rep(NA_character_, length(taxids))
  names(fam) <- taxids
  
  failed_taxids <- character(0)
  chunks <- split(taxids, ceiling(seq_along(taxids) / chunk_size))
  
  parse_classification_result <- function(cl_obj) {
    if (is.null(cl_obj)) return(NULL)
    
    for (tid in names(cl_obj)) {
      entry <- cl_obj[[tid]]
      
      if (is.data.frame(entry)) {
        fam_row <- entry[entry$rank == "family", , drop = FALSE]
        
        if (nrow(fam_row) > 0) {
          fam[tid] <<- fam_row$name[1]
        } else {
          fam[tid] <<- NA_character_
        }
      } else {
        fam[tid] <<- NA_character_
      }
    }
  }
  
  for (i in seq_along(chunks)) {
    this_chunk <- chunks[[i]]
    
    message(
      "Resolving chunk ", i, " of ", length(chunks),
      " (", length(this_chunk), " taxids)"
    )
    
    cl <- tryCatch(
      taxize::classification(this_chunk, db = "ncbi"),
      error = function(e) {
        message("Chunk ", i, " failed: ", conditionMessage(e))
        failed_taxids <<- c(failed_taxids, this_chunk)
        return(NULL)
      }
    )
    
    parse_classification_result(cl)
    Sys.sleep(sleep_sec)
  }
  
  failed_taxids <- unique(failed_taxids)
  
  if (retry_singletons && length(failed_taxids) > 0) {
    message("Retrying ", length(failed_taxids), " failed taxids one at a time.")
    
    still_failed <- character(0)
    
    for (j in seq_along(failed_taxids)) {
      tid <- failed_taxids[j]
      
      message("Retrying failed taxid ", j, " of ", length(failed_taxids), ": ", tid)
      
      cl_single <- tryCatch(
        taxize::classification(tid, db = "ncbi"),
        error = function(e) {
          message("Taxid ", tid, " failed: ", conditionMessage(e))
          still_failed <<- c(still_failed, tid)
          return(NULL)
        }
      )
      
      parse_classification_result(cl_single)
      Sys.sleep(sleep_sec)
    }
    
    failed_taxids <- unique(still_failed)
  }
  
  attr(fam, "failed_taxids") <- failed_taxids
  fam
}

family_lookup <- get_family_safe(
  taxid_vec,
  chunk_size = 50,
  sleep_sec = 0.4,
  retry_singletons = TRUE
)

failed_taxids <- attr(family_lookup, "failed_taxids")

if (length(failed_taxids) > 0) {
  readr::write_tsv(
    tibble::tibble(taxid = failed_taxids),
    file.path(family_dirs$results, "failed_family_lookup_taxids.tsv")
  )
  
  message(
    "Some taxids failed. Failed taxids written to: ",
    file.path(family_dirs$results, "failed_family_lookup_taxids.tsv")
  )
}

family_lookup_tbl <- tibble::tibble(
  taxid = names(family_lookup),
  family = unname(family_lookup)
)

readr::write_tsv(
  family_lookup_tbl,
  file.path(family_dirs$data, "taxid_to_family_lookup.tsv")
)

cat(
  "Resolved families for",
  sum(!is.na(family_lookup_tbl$family)),
  "of",
  nrow(family_lookup_tbl),
  "taxids.\n"
)

long_tax <- long_tax |>
  dplyr::mutate(taxid = as.character(taxid)) |>
  dplyr::left_join(family_lookup_tbl, by = "taxid")

long_tax$family[is.na(long_tax$family)] <- "unknown_family"

readr::write_tsv(
  long_tax,
  file.path(family_dirs$data, "abundance_long_with_taxid_and_family.tsv")
)

family_sample <- long_tax |>
  dplyr::mutate(status = ifelse(endo, "endogenous", "free")) |>
  dplyr::group_by(sample, family, status) |>
  dplyr::summarise(reads = sum(clade_reads), .groups = "drop")

family_totals <- family_sample |>
  dplyr::group_by(family, status) |>
  dplyr::summarise(total_reads = sum(reads), .groups = "drop") |>
  tidyr::pivot_wider(
    names_from = status,
    values_from = total_reads,
    values_fill = 0
  )

if (!"endogenous" %in% colnames(family_totals)) family_totals[["endogenous"]] <- 0
if (!"free" %in% colnames(family_totals)) family_totals[["free"]] <- 0

family_totals[["endo_reads"]] <- as.numeric(family_totals[["endogenous"]])
family_totals[["nonendo_reads"]] <- as.numeric(family_totals[["free"]])

family_totals <- family_totals |>
  dplyr::select(-dplyr::any_of(c("endogenous", "free"))) |>
  dplyr::mutate(
    total_reads = endo_reads + nonendo_reads,
    prop_endo = ifelse(total_reads > 0, endo_reads / total_reads, 0)
  ) |>
  dplyr::arrange(desc(total_reads))

readr::write_tsv(family_totals, file.path(family_dirs$data, "family_endogenization_totals.tsv"))

global_endo <- sum(family_totals$endo_reads, na.rm = TRUE)
global_nonendo <- sum(family_totals$nonendo_reads, na.rm = TRUE)

stats_list <- purrr::map_dfr(seq_len(nrow(family_totals)), function(i) {
  row <- family_totals[i, ]
  
  a <- row$endo_reads
  b <- row$nonendo_reads
  c <- global_endo - a
  d <- global_nonendo - b
  
  tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  
  ft <- tryCatch(
    fisher.test(tab),
    error = function(e) {
      fisher.test(tab + 0.5)
    }
  )
  
  tibble::tibble(
    family = row$family,
    odds_ratio = as.numeric(ft$estimate),
    p_value = ft$p.value,
    endo_reads = a,
    nonendo_reads = b,
    total_reads = a + b
  )
})

stats_list <- stats_list |>
  dplyr::mutate(padj = p.adjust(p_value, method = "BH")) |>
  dplyr::arrange(padj)

readr::write_tsv(stats_list, file.path(family_dirs$results, "family_endogenization_fisher_results.tsv"))

summary_out <- stats_list |>
  dplyr::left_join(family_totals |> dplyr::select(family, prop_endo), by = "family") |>
  dplyr::mutate(significant = padj < 0.05)

readr::write_tsv(
  summary_out,
  file.path(family_dirs$results, "family_endogenization_summary_with_significance.tsv")
)
# ----------------------------
# HOST-FACETED FAMILY-LEVEL ENDOGENIZATION COMPARISON
# ----------------------------

cat("Running host-faceted family-level endogenization analysis...\n")

host_family_dirs <- list(
  base = file.path(family_dirs$base, "host_faceted"),
  data = file.path(family_dirs$base, "host_faceted", "data"),
  results = file.path(family_dirs$base, "host_faceted", "results"),
  figures = file.path(family_dirs$base, "host_faceted", "figures")
)

purrr::walk(host_family_dirs, ~dir.create(.x, showWarnings = FALSE, recursive = TRUE))

# Settings
top_n_families <- 15
other_label <- "Other"
min_reads_per_host_family <- 10
include_other_in_forest <- TRUE

status_colors <- c(
  "Endogenous" = "#de2d26",
  "Free virus" = "#4575b4"
)

host_order <- c(
  "Aedes albopictus",
  "Culex pipiens",
  "Culiseta melanura"
)

# This section uses objects already created earlier in the script:
# long_tax: sample, taxid, taxon, clade_reads, endo, family
# meta_for_plot: sample, host

required_long_tax_cols <- c("sample", "family", "clade_reads", "endo")
missing_long_tax_cols <- setdiff(required_long_tax_cols, colnames(long_tax))

if (length(missing_long_tax_cols) > 0) {
  stop(
    "long_tax is missing required columns for host-faceted family analysis: ",
    paste(missing_long_tax_cols, collapse = ", ")
  )
}

required_meta_cols <- c("sample", "host")
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta_for_plot))

if (length(missing_meta_cols) > 0) {
  stop(
    "meta_for_plot is missing required columns for host-faceted family analysis: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

# Join host metadata to family-level viral read calls.
# Unknown/unresolved families are collapsed into "Other" here so they do not appear
# as a separate unknown_family category in either host-faceted plot.
host_family_long <- long_tax |>
  dplyr::left_join(
    meta_for_plot |> dplyr::select(sample, host),
    by = "sample"
  ) |>
  dplyr::mutate(
    host = factor(host, levels = host_order),
    family = dplyr::case_when(
      is.na(family) ~ other_label,
      family == "" ~ other_label,
      family %in% c(
        "unknown_family",
        "Unknown family",
        "unknown",
        "Unknown",
        "Unassigned family",
        "unassigned_family",
        "NA"
      ) ~ other_label,
      TRUE ~ as.character(family)
    ),
    status = dplyr::if_else(endo, "Endogenous", "Free virus"),
    status = factor(status, levels = c("Free virus", "Endogenous")),
    clade_reads = as.numeric(clade_reads)
  ) |>
  dplyr::filter(
    !is.na(host),
    !is.na(family),
    !is.na(status),
    !is.na(clade_reads),
    clade_reads > 0
  )

readr::write_tsv(
  host_family_long,
  file.path(host_family_dirs$data, "host_family_long_input.tsv")
)

# Collapse to family x host x status counts
host_family_counts_raw <- host_family_long |>
  dplyr::group_by(host, family, status) |>
  dplyr::summarise(
    reads = sum(clade_reads, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_tsv(
  host_family_counts_raw,
  file.path(host_family_dirs$data, "host_family_counts_raw.tsv")
)

# Select top families globally so all host facets use the same family set/order.
# "Other" is excluded from top-family selection and then added back after collapse.
top_families <- host_family_counts_raw |>
  dplyr::filter(family != other_label) |>
  dplyr::group_by(family) |>
  dplyr::summarise(
    total_reads = sum(reads, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(total_reads)) |>
  dplyr::slice_head(n = top_n_families) |>
  dplyr::pull(family)

host_family_counts <- host_family_counts_raw |>
  dplyr::mutate(
    family_display = dplyr::if_else(family %in% top_families, family, other_label)
  ) |>
  dplyr::group_by(host, family_display, status) |>
  dplyr::summarise(
    reads = sum(reads, na.rm = TRUE),
    .groups = "drop"
  )

family_order <- host_family_counts |>
  dplyr::group_by(family_display) |>
  dplyr::summarise(
    total_reads = sum(reads, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(total_reads)) |>
  dplyr::pull(family_display)

if (other_label %in% family_order) {
  family_order <- c(setdiff(family_order, other_label), other_label)
}

# Expand missing combinations so all host facets show the same family/status layout
host_family_counts_complete <- tidyr::expand_grid(
  host = factor(host_order, levels = host_order),
  family_display = family_order,
  status = factor(c("Free virus", "Endogenous"), levels = c("Free virus", "Endogenous"))
) |>
  dplyr::left_join(
    host_family_counts,
    by = c("host", "family_display", "status")
  ) |>
  dplyr::mutate(
    reads = tidyr::replace_na(reads, 0),
    family_display = factor(family_display, levels = family_order),
    status = factor(status, levels = c("Free virus", "Endogenous"))
  )

# Summary proportions for stacked bar plot
host_family_props <- host_family_counts_complete |>
  dplyr::group_by(host, family_display) |>
  dplyr::mutate(
    total_reads_family_host = sum(reads, na.rm = TRUE),
    proportion = dplyr::if_else(
      total_reads_family_host > 0,
      reads / total_reads_family_host,
      0
    )
  ) |>
  dplyr::ungroup()

readr::write_tsv(
  host_family_props,
  file.path(host_family_dirs$results, "host_faceted_family_proportions.tsv")
)

# ----------------------------
# HOST-SPECIFIC FISHER TESTS
# ----------------------------

sig_label <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p <= 0.0001 ~ "****",
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

compute_host_family_fisher <- function(a, b, c, d) {
  tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  
  ft <- tryCatch(
    stats::fisher.test(tab),
    error = function(e) NULL
  )
  
  fisher_p <- if (!is.null(ft)) ft$p.value else NA_real_
  fisher_or <- if (!is.null(ft)) as.numeric(ft$estimate) else NA_real_
  
  # Haldane-Anscombe correction for stable log2 OR and CI plotting
  a2 <- a + 0.5
  b2 <- b + 0.5
  c2 <- c + 0.5
  d2 <- d + 0.5
  
  odds_ratio_corrected <- (a2 * d2) / (b2 * c2)
  log_or <- log(odds_ratio_corrected)
  se_log <- sqrt((1 / a2) + (1 / b2) + (1 / c2) + (1 / d2))
  
  tibble::tibble(
    fisher_p_value = fisher_p,
    fisher_odds_ratio = fisher_or,
    odds_ratio_corrected = odds_ratio_corrected,
    log2_odds_ratio = log_or / log(2),
    conf_low_log2 = (log_or - 1.96 * se_log) / log(2),
    conf_high_log2 = (log_or + 1.96 * se_log) / log(2)
  )
}

host_family_wide <- host_family_counts_complete |>
  tidyr::pivot_wider(
    names_from = status,
    values_from = reads,
    values_fill = 0
  ) |>
  dplyr::rename(
    free_reads = `Free virus`,
    endo_reads = Endogenous
  ) |>
  dplyr::mutate(
    total_reads = endo_reads + free_reads
  )

host_totals <- host_family_wide |>
  dplyr::group_by(host) |>
  dplyr::summarise(
    host_endo_reads = sum(endo_reads, na.rm = TRUE),
    host_free_reads = sum(free_reads, na.rm = TRUE),
    .groups = "drop"
  )

host_fisher_input <- host_family_wide |>
  dplyr::left_join(host_totals, by = "host") |>
  dplyr::mutate(
    other_endo_reads = host_endo_reads - endo_reads,
    other_free_reads = host_free_reads - free_reads
  ) |>
  dplyr::filter(total_reads >= min_reads_per_host_family)

if (!include_other_in_forest) {
  host_fisher_input <- host_fisher_input |>
    dplyr::filter(family_display != other_label)
}

host_fisher_results <- host_fisher_input |>
  dplyr::rowwise() |>
  dplyr::mutate(
    fisher_stats = list(
      compute_host_family_fisher(
        a = endo_reads,
        b = free_reads,
        c = other_endo_reads,
        d = other_free_reads
      )
    )
  ) |>
  tidyr::unnest(fisher_stats) |>
  dplyr::ungroup() |>
  dplyr::group_by(host) |>
  dplyr::mutate(
    p_adj_BH_within_host = p.adjust(fisher_p_value, method = "BH"),
    significance_within_host = sig_label(p_adj_BH_within_host)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    p_adj_BH_global = p.adjust(fisher_p_value, method = "BH"),
    significance_global = sig_label(p_adj_BH_global),
    host = factor(host, levels = host_order),
    family_display = factor(family_display, levels = rev(family_order))
  ) |>
  dplyr::arrange(host, p_adj_BH_within_host)

readr::write_tsv(
  host_fisher_results,
  file.path(host_family_dirs$results, "host_faceted_family_fisher_results.tsv")
)

# Summary table
host_family_summary <- host_family_wide |>
  dplyr::mutate(
    prop_endo = dplyr::if_else(total_reads > 0, endo_reads / total_reads, NA_real_),
    prop_free = dplyr::if_else(total_reads > 0, free_reads / total_reads, NA_real_)
  ) |>
  dplyr::left_join(
    host_fisher_results |>
      dplyr::select(
        host,
        family_display,
        log2_odds_ratio,
        conf_low_log2,
        conf_high_log2,
        fisher_p_value,
        p_adj_BH_within_host,
        significance_within_host,
        p_adj_BH_global,
        significance_global
      ),
    by = c("host", "family_display")
  ) |>
  dplyr::arrange(host, dplyr::desc(total_reads))

readr::write_tsv(
  host_family_summary,
  file.path(host_family_dirs$results, "host_faceted_family_summary.tsv")
)

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  
  openxlsx::addWorksheet(wb, "host_family_summary")
  openxlsx::writeData(wb, "host_family_summary", host_family_summary)
  
  openxlsx::addWorksheet(wb, "fisher_results")
  openxlsx::writeData(wb, "fisher_results", host_fisher_results)
  
  openxlsx::addWorksheet(wb, "proportions")
  openxlsx::writeData(wb, "proportions", host_family_props)
  
  openxlsx::saveWorkbook(
    wb,
    file.path(host_family_dirs$results, "host_faceted_family_endogenization_summary.xlsx"),
    overwrite = TRUE
  )
}

# ----------------------------
# PLOT 1: HOST-FACETED STACKED BAR PLOT
# ----------------------------

p_host_bar <- ggplot(
  host_family_props,
  aes(x = family_display, y = proportion, fill = status)
) +
  geom_col(width = 0.9) +
  facet_wrap(~ host, nrow = 1) +
  scale_fill_manual(values = status_colors, drop = FALSE) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Proportion endogenous vs free viruses by host",
    x = "Viral family",
    y = "Proportion",
    fill = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92", color = "black"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(
  file.path(host_family_dirs$figures, "host_faceted_endogenous_vs_free_proportions.pdf"),
  p_host_bar,
  width = 16,
  height = 6.5,
  bg = "white"
)

ggsave(
  file.path(host_family_dirs$figures, "host_faceted_endogenous_vs_free_proportions.png"),
  p_host_bar,
  width = 16,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    file.path(host_family_dirs$figures, "host_faceted_endogenous_vs_free_proportions.svg"),
    p_host_bar,
    width = 16,
    height = 6.5,
    device = svglite::svglite,
    bg = "white"
  )
}

# ----------------------------
# PLOT 2: HOST-FACETED FOREST PLOT
# ----------------------------

forest_plot_df <- host_fisher_results |>
  dplyr::filter(
    !is.na(log2_odds_ratio),
    !is.na(conf_low_log2),
    !is.na(conf_high_log2)
  )

if (nrow(forest_plot_df) > 0) {
  
  forest_range <- range(
    c(forest_plot_df$conf_low_log2, forest_plot_df$conf_high_log2),
    na.rm = TRUE
  )
  
  forest_pad <- diff(forest_range) * 0.08
  if (!is.finite(forest_pad) || forest_pad == 0) {
    forest_pad <- 1
  }
  
  forest_limits <- c(forest_range[1] - forest_pad, forest_range[2] + forest_pad)
  
  p_host_forest <- ggplot(
    forest_plot_df,
    aes(x = log2_odds_ratio, y = family_display)
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey45",
      linewidth = 0.7
    ) +
    geom_segment(
      aes(
        x = conf_low_log2,
        xend = conf_high_log2,
        yend = family_display
      ),
      linewidth = 0.6,
      color = "black"
    ) +
    geom_point(
      aes(size = total_reads),
      shape = 21,
      fill = "white",
      color = "black",
      stroke = 0.9
    ) +
    facet_wrap(~ host, nrow = 1, scales = "fixed") +
    scale_size_continuous(
      range = c(2.2, 9),
      labels = scales::comma
    ) +
    coord_cartesian(xlim = forest_limits, clip = "off") +
    labs(
      title = "Endogenization enrichment by viral family within host",
      x = "log2 odds ratio (endogenous vs free virus)",
      y = "Viral family",
      size = "Total reads"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92", color = "black"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  
  ggsave(
    file.path(host_family_dirs$figures, "host_faceted_endogenization_forest_plot.pdf"),
    p_host_forest,
    width = 16,
    height = 8,
    bg = "white"
  )
  
  ggsave(
    file.path(host_family_dirs$figures, "host_faceted_endogenization_forest_plot.png"),
    p_host_forest,
    width = 16,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(
      file.path(host_family_dirs$figures, "host_faceted_endogenization_forest_plot.svg"),
      p_host_forest,
      width = 16,
      height = 8,
      device = svglite::svglite,
      bg = "white"
    )
  }
  
} else {
  warning("No valid host-specific Fisher results were available for forest plotting.")
}

cat("\nHost-faceted family-level endogenization analysis complete.\n")
cat("Outputs written to:\n")
cat("  ", host_family_dirs$data, "\n")
cat("  ", host_family_dirs$results, "\n")
cat("  ", host_family_dirs$figures, "\n")
# ----------------------------
# FAMILY BAR PLOTS
# ----------------------------

plot_df <- family_totals |>
  dplyr::select(family, endo_reads, nonendo_reads) |>
  tidyr::pivot_longer(
    cols = c(endo_reads, nonendo_reads),
    names_to = "status",
    values_to = "reads"
  ) |>
  dplyr::mutate(
    status = dplyr::recode(
      status,
      endo_reads = "Endogenous",
      nonendo_reads = "Free virus"
    ),
    family = factor(family, levels = family_totals$family)
  )

p_abs <- ggplot(plot_df, aes(x = family, y = reads, fill = status)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Virus family abundance: endogenous vs free",
    x = "Viral family",
    y = "Total reads"
  ) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette = "Set1")

ggsave(
  file.path(family_dirs$figures, "family_endogenization_absolute_barplot.pdf"),
  p_abs,
  width = 10,
  height = 6,
  bg = "white"
)

p_prop <- ggplot(plot_df, aes(x = family, y = reads, fill = status)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Proportion endogenous vs free by viral family",
    x = "Viral family",
    y = "Proportion"
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette = "Set1")

ggsave(
  file.path(family_dirs$figures, "family_endogenization_proportion_barplot.pdf"),
  p_prop,
  width = 10,
  height = 6,
  bg = "white"
)

# ----------------------------
# FAMILY FOREST PLOT HELPERS
# ----------------------------

compute_or_ci <- function(a, b, c, d) {
  tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  
  ft <- tryCatch(fisher.test(tab), error = function(e) NULL)
  if (is.null(ft) || is.null(ft$estimate) || is.na(ft$estimate)) {
    ft <- tryCatch(fisher.test(tab + 0.5), error = function(e) NULL)
  }
  
  if (is.null(ft)) {
    return(tibble::tibble(
      odds_ratio = NA_real_,
      log2_or = NA_real_,
      log2_ci_low = NA_real_,
      log2_ci_high = NA_real_,
      p_value = NA_real_
    ))
  }
  
  or <- as.numeric(ft$estimate)
  ci <- as.numeric(ft$conf.int)
  p <- as.numeric(ft$p.value)
  
  tibble::tibble(
    odds_ratio = or,
    log2_or = ifelse(is.finite(or) && or > 0, log2(or), NA_real_),
    log2_ci_low = ifelse(is.finite(ci[1]) && ci[1] > 0, log2(ci[1]), NA_real_),
    log2_ci_high = ifelse(is.finite(ci[2]) && ci[2] > 0, log2(ci[2]), NA_real_),
    p_value = p
  )
}

forest_all <- family_totals |>
  dplyr::rowwise() |>
  dplyr::mutate(
    other_endo = global_endo - endo_reads,
    other_nonendo = global_nonendo - nonendo_reads,
    ci = list(compute_or_ci(endo_reads, nonendo_reads, other_endo, other_nonendo)),
    odds_ratio = ci$odds_ratio,
    log2_or = ci$log2_or,
    log2_ci_low = ci$log2_ci_low,
    log2_ci_high = ci$log2_ci_high,
    p_value = ci$p_value
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(padj = p.adjust(p_value, method = "BH"))

readr::write_tsv(forest_all, file.path(family_dirs$data, "family_endogenization_forest_data.tsv"))

# ----------------------------
# COLLAPSE LOW-COUNT FAMILIES INTO "Other"
# ----------------------------

top_n <- 15

family_totals_collapsed_base <- family_totals |>
  dplyr::mutate(
    family = dplyr::case_when(
      is.na(family) ~ "Other",
      family %in% c("unknown_family", "Unassigned family") ~ "Other",
      TRUE ~ family
    )
  ) |>
  dplyr::group_by(family) |>
  dplyr::summarise(
    endo_reads = sum(endo_reads, na.rm = TRUE),
    nonendo_reads = sum(nonendo_reads, na.rm = TRUE),
    total_reads = sum(total_reads, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(prop_endo = ifelse(total_reads > 0, endo_reads / total_reads, 0))

global_endo_collapsed <- sum(family_totals_collapsed_base$endo_reads, na.rm = TRUE)
global_nonendo_collapsed <- sum(family_totals_collapsed_base$nonendo_reads, na.rm = TRUE)
global_total <- global_endo_collapsed + global_nonendo_collapsed
global_prop_endo <- ifelse(global_total > 0, global_endo_collapsed / global_total, NA_real_)

top_fams <- family_totals_collapsed_base |>
  dplyr::filter(family != "Other") |>
  dplyr::arrange(desc(total_reads)) |>
  dplyr::slice_head(n = top_n) |>
  dplyr::pull(family)

collapsed <- family_totals_collapsed_base |>
  dplyr::mutate(family_group = ifelse(family %in% top_fams, as.character(family), "Other")) |>
  dplyr::group_by(family_group) |>
  dplyr::summarise(
    endo_reads = sum(endo_reads, na.rm = TRUE),
    nonendo_reads = sum(nonendo_reads, na.rm = TRUE),
    total_reads = sum(total_reads, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    other_endo = pmax(0, global_endo_collapsed - endo_reads),
    other_nonendo = pmax(0, global_nonendo_collapsed - nonendo_reads),
    ci = list(compute_or_ci(endo_reads, nonendo_reads, other_endo, other_nonendo)),
    odds_ratio = ci$odds_ratio,
    log2_or = ci$log2_or,
    log2_ci_low = ci$log2_ci_low,
    log2_ci_high = ci$log2_ci_high,
    p_value = ci$p_value
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    padj = p.adjust(p_value, method = "BH"),
    family = family_group,
    prop_endo = ifelse(total_reads > 0, endo_reads / total_reads, 0),
    expected_endo_reads = total_reads * global_prop_endo,
    observed_expected_ratio = ifelse(expected_endo_reads > 0, endo_reads / expected_endo_reads, NA_real_),
    log2_obs_exp = ifelse(
      is.finite(observed_expected_ratio) & observed_expected_ratio > 0,
      log2(observed_expected_ratio),
      NA_real_
    )
  ) |>
  dplyr::select(
    family,
    endo_reads,
    nonendo_reads,
    total_reads,
    prop_endo,
    odds_ratio,
    log2_or,
    log2_ci_low,
    log2_ci_high,
    p_value,
    padj,
    expected_endo_reads,
    observed_expected_ratio,
    log2_obs_exp
  ) |>
  dplyr::arrange(log2_or)

readr::write_tsv(
  collapsed,
  file.path(family_dirs$data, paste0("family_forest_collapsed_top", top_n, "_data.tsv"))
)

# ----------------------------
# COLLAPSED FOREST PLOT
# ----------------------------

collapsed_plot_df <- collapsed |>
  dplyr::filter(
    is.finite(log2_or),
    is.finite(log2_ci_low),
    is.finite(log2_ci_high)
  ) |>
  dplyr::arrange(log2_or) |>
  dplyr::mutate(family = factor(family, levels = family))

p_forest_collapsed <- ggplot(collapsed_plot_df, aes(x = log2_or, y = family)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.7
  ) +
  geom_errorbarh(
    aes(xmin = log2_ci_low, xmax = log2_ci_high),
    height = 0.22,
    linewidth = 0.9,
    color = "black"
  ) +
  geom_point(
    aes(size = total_reads),
    shape = 21,
    fill = NA,
    color = "black",
    stroke = 1.1
  ) +
  scale_size_continuous(
    range = c(3, 9),
    breaks = scales::pretty_breaks(n = 4),
    labels = scales::label_comma(),
    name = "Total reads"
  ) +
  scale_x_continuous(
    name = "log2 odds ratio (endogenous vs free virus)",
    breaks = scales::pretty_breaks(n = 7)
  ) +
  labs(
    title = "Endogenization enrichment by viral family",
    y = "Viral family"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray82"),
    panel.grid.major.x = element_line(color = "gray85"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  file.path(family_dirs$figures, paste0("family_forest_collapsed_top", top_n, ".pdf")),
  p_forest_collapsed,
  width = 10.5,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(family_dirs$figures, paste0("family_forest_collapsed_top", top_n, ".png")),
  p_forest_collapsed,
  width = 10.5,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

# ----------------------------
# COLLAPSED STACKED PROPORTION BARPLOT
# ----------------------------

collapsed_bar_df <- collapsed |>
  dplyr::arrange(prop_endo) |>
  dplyr::mutate(family = factor(family, levels = family)) |>
  dplyr::select(family, endo_reads, nonendo_reads, total_reads) |>
  tidyr::pivot_longer(
    cols = c(endo_reads, nonendo_reads),
    names_to = "status",
    values_to = "reads"
  ) |>
  dplyr::mutate(
    status = dplyr::recode(
      status,
      endo_reads = "Endogenous",
      nonendo_reads = "Free virus"
    )
  )

p_collapsed_bar <- ggplot(collapsed_bar_df, aes(x = family, y = reads, fill = status)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(
    values = c(
      "Endogenous" = "#D73027",
      "Free virus" = "#4575B4"
    ),
    name = NULL
  ) +
  labs(
    title = paste0("Proportion endogenous vs free viruses (Top ", top_n, " + Other)"),
    x = "Viral family",
    y = "Proportion"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  file.path(family_dirs$figures, paste0("family_collapsed_proportion_top", top_n, ".pdf")),
  p_collapsed_bar,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(family_dirs$figures, paste0("family_collapsed_proportion_top", top_n, ".png")),
  p_collapsed_bar,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)

# ----------------------------
# OBSERVED VS EXPECTED ENDO ENRICHMENT
# ----------------------------

obs_exp_df <- collapsed |>
  dplyr::filter(is.finite(log2_obs_exp)) |>
  dplyr::arrange(log2_obs_exp) |>
  dplyr::mutate(family = factor(family, levels = family))

p_obs_exp <- ggplot(obs_exp_df, aes(x = log2_obs_exp, y = family)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray35", linewidth = 0.7) +
  geom_point(
    aes(size = total_reads),
    shape = 21,
    fill = NA,
    color = "black",
    stroke = 1.1
  ) +
  scale_size_continuous(
    range = c(3, 9),
    breaks = scales::pretty_breaks(n = 4),
    labels = scales::label_comma(),
    name = "Total reads"
  ) +
  scale_x_continuous(
    name = "log2 observed / expected endogenous reads",
    breaks = scales::pretty_breaks(n = 7)
  ) +
  labs(
    title = "Observed vs expected endogenous enrichment by viral family",
    y = "Viral family"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray82"),
    panel.grid.major.x = element_line(color = "gray85"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  file.path(family_dirs$figures, paste0("family_observed_expected_top", top_n, ".pdf")),
  p_obs_exp,
  width = 10.5,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(family_dirs$figures, paste0("family_observed_expected_top", top_n, ".png")),
  p_obs_exp,
  width = 10.5,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

# ----------------------------
# OUTPUT SUMMARY
# ----------------------------

writeLines(
  c(
    paste0("forest_collapsed: ", file.path(family_dirs$figures, paste0("family_forest_collapsed_top", top_n, ".pdf"))),
    paste0("collapsed_bar: ", file.path(family_dirs$figures, paste0("family_collapsed_proportion_top", top_n, ".pdf"))),
    paste0("observed_expected: ", file.path(family_dirs$figures, paste0("family_observed_expected_top", top_n, ".pdf")))
  ),
  con = file.path(family_dirs$results, "family_forest_and_collapsed_outputs.txt")
)

# ----------------------------
# COVERAGE ANALYSIS
# ----------------------------

coverage_outdir <- coverage_dirs$base
coverage_data_dir <- coverage_dirs$data
coverage_results_dir <- coverage_dirs$results
coverage_figures_dir <- coverage_dirs$figures

df <- readr::read_csv(summary_file, show_col_types = FALSE)

# Make column names syntactically safe while preserving their meaning.
colnames(df) <- make.names(colnames(df))

cat("\nColumns detected in ", summary_file, " after make.names():\n", sep = "")
print(colnames(df))

# Expected columns in summary_complete.csv after make.names():
# Sample
# Mosquito_reads
# Percent_mosquito
# Non.mosquito_reads
# Total_reads
# Viral_reads
# Percent_viral
# Average_coverage
# Normalized.coverage
# Enriched
# If Host is absent from this file, host will be inferred from metadata_file when possible.

find_first_existing_col <- function(df, candidates, required_name) {
  hit <- candidates[candidates %in% colnames(df)]
  
  if (length(hit) == 0) {
    stop(
      "Could not find required column for ", required_name, ".\n",
      "Tried: ", paste(candidates, collapse = ", "), "\n",
      "Available columns are:\n",
      paste(colnames(df), collapse = ", ")
    )
  }
  
  hit[1]
}

sample_col <- find_first_existing_col(
  df,
  c("Sample", "sample"),
  "sample"
)

host_col <- c("Host", "host")[c("Host", "host") %in% colnames(df)]
host_col <- if (length(host_col) > 0) host_col[1] else NA_character_

percent_mosquito_col <- find_first_existing_col(
  df,
  c("Percent_mosquito", "Percent.mosquito", "percent_mosquito"),
  "percent mosquito"
)

viral_reads_col <- find_first_existing_col(
  df,
  c("Viral_reads", "Viral.reads", "viral_reads"),
  "viral reads"
)

other_reads_col <- find_first_existing_col(
  df,
  c("Non.mosquito_reads", "Non_mosquito_reads", "Other_reads", "other_reads"),
  "non-mosquito/other reads"
)

average_coverage_col <- find_first_existing_col(
  df,
  c("Average_coverage", "Average.coverage", "average_coverage"),
  "average mosquito genome coverage"
)

normalized_coverage_col <- find_first_existing_col(
  df,
  c(
    "Mosquito.genome.coverage.per.million.mosquito.reads",
    "Mosquito_genome_coverage_per_million_mosquito_reads",
    "Normalized.coverage",
    "Normalized_coverage",
    "norm_cov"
  ),
  "mosquito genome coverage per million mosquito reads"
)

enrichment_col <- find_first_existing_col(
  df,
  c("Enriched", "enriched", "Enrichment", "enrichment"),
  "enrichment"
)

# Standardize names used in the plotting/statistical code.
if (!is.na(host_col)) {
  df <- df |>
    dplyr::rename(
      sample = dplyr::all_of(sample_col),
      host = dplyr::all_of(host_col),
      percent_mosquito = dplyr::all_of(percent_mosquito_col),
      viral_reads = dplyr::all_of(viral_reads_col),
      other_reads = dplyr::all_of(other_reads_col),
      average_coverage = dplyr::all_of(average_coverage_col),
      normalized_mosquito_coverage = dplyr::all_of(normalized_coverage_col),
      enrichment = dplyr::all_of(enrichment_col)
    )
} else {
  df <- df |>
    dplyr::rename(
      sample = dplyr::all_of(sample_col),
      percent_mosquito = dplyr::all_of(percent_mosquito_col),
      viral_reads = dplyr::all_of(viral_reads_col),
      other_reads = dplyr::all_of(other_reads_col),
      average_coverage = dplyr::all_of(average_coverage_col),
      normalized_mosquito_coverage = dplyr::all_of(normalized_coverage_col),
      enrichment = dplyr::all_of(enrichment_col)
    ) |>
    dplyr::left_join(meta_for_plot, by = "sample")
}

if (!"host" %in% colnames(df)) {
  df$host <- NA_character_
}

# Coerce enrichment robustly in case it is read as logical, character, or numeric.
df <- df |>
  dplyr::mutate(
    enrichment = dplyr::case_when(
      enrichment %in% c(TRUE, "TRUE", "True", "true", "T", "t", "1", 1) ~ TRUE,
      enrichment %in% c(FALSE, "FALSE", "False", "false", "F", "f", "0", 0) ~ FALSE,
      TRUE ~ NA
    ),
    enrichment = factor(enrichment, levels = c(FALSE, TRUE), labels = c("Unenriched", "Enriched")),
    viral_log10 = log10(viral_reads + 1),
    other_log10 = log10(other_reads + 1)
  )

if (any(is.na(df$enrichment))) {
  warning("Some rows have missing/unclear Enriched values after coercion.")
}

coverage_plot_df <- df |>
  dplyr::select(
    sample,
    host,
    enrichment,
    average_coverage,
    normalized_mosquito_coverage,
    percent_mosquito,
    viral_log10,
    other_log10
  ) |>
  tidyr::pivot_longer(
    cols = c(
      average_coverage,
      normalized_mosquito_coverage,
      percent_mosquito,
      viral_log10,
      other_log10
    ),
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::mutate(
    metric = dplyr::recode(
      metric,
      average_coverage = "Average mosquito genome coverage",
      normalized_mosquito_coverage = "Mosquito genome coverage per million mosquito reads",
      percent_mosquito = "Percent mosquito",
      viral_log10 = "Viral reads (log10 + 1)",
      other_log10 = "Non-mosquito reads (log10 + 1)"
    ),
    metric = factor(
      metric,
      levels = c(
        "Average mosquito genome coverage",
        "Mosquito genome coverage per million mosquito reads",
        "Percent mosquito",
        "Viral reads (log10 + 1)",
        "Non-mosquito reads (log10 + 1)"
      )
    )
  )

readr::write_csv(
  coverage_plot_df,
  file.path(coverage_data_dir, "coverage_plot_long.csv")
)

# Wilcoxon tests comparing enriched vs unenriched samples for each metric.
coverage_stats_df <- coverage_plot_df |>
  dplyr::filter(!is.na(value), !is.na(enrichment)) |>
  dplyr::group_by(metric) |>
  rstatix::wilcox_test(value ~ enrichment) |>
  rstatix::adjust_pvalue(method = "BH") |>
  rstatix::add_significance("p.adj") |>
  dplyr::left_join(
    coverage_plot_df |>
      dplyr::filter(!is.na(value), !is.na(enrichment)) |>
      dplyr::group_by(metric) |>
      rstatix::wilcox_effsize(value ~ enrichment, ci = FALSE),
    by = "metric"
  ) |>
  dplyr::left_join(
    coverage_plot_df |>
      dplyr::group_by(metric) |>
      dplyr::summarise(
        ymax = max(value, na.rm = TRUE),
        ymin = min(value, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "metric"
  ) |>
  dplyr::mutate(
    y.label = ymax + (ymax - ymin) * 0.08,
    label_text = paste0(p.adj.signif, "   r = ", sprintf("%.2f", effsize)),
    x = 1.5
  )

readr::write_csv(
  coverage_stats_df,
  file.path(coverage_results_dir, "coverage_analysis_wilcox_results.csv")
)

coverage_stats_summary <- coverage_stats_df |>
  dplyr::select(dplyr::any_of(c(
    "metric",
    ".y.",
    "group1",
    "group2",
    "n1",
    "n2",
    "statistic",
    "p",
    "p.adj",
    "p.adj.signif",
    "effsize"
  ))) |>
  dplyr::rename(
    wilcox_statistic = statistic,
    raw_p_value = p,
    adjusted_p_value = p.adj,
    significance = p.adj.signif,
    rank_biserial_r = effsize
  )

readr::write_csv(
  coverage_stats_summary,
  file.path(coverage_results_dir, "coverage_analysis_wilcox_summary.csv")
)

p_coverage <- ggplot(coverage_plot_df, aes(x = enrichment, y = value, color = enrichment)) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.85) +
  geom_boxplot(alpha = 0.25, outlier.shape = NA) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  geom_text(
    data = coverage_stats_df,
    aes(x = x, y = y.label, label = label_text),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "Effect of enrichment on mosquito coverage and read composition",
    x = "Sample type",
    y = "Measurement"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(coverage_figures_dir, "enrichment_effects_panels_with_effectsize.pdf"),
  p_coverage,
  width = 11,
  height = 9,
  bg = "white"
)

ggsave(
  file.path(coverage_figures_dir, "enrichment_effects_panels_with_effectsize.png"),
  p_coverage,
  width = 11,
  height = 9,
  dpi = 300,
  bg = "white"
)
# ----------------------------
# TEST CODE: ADONIS WITH MARGINAL EFFECTS
# ----------------------------

run_adonis_margin <- function(dist_obj, meta_tbl) {
  df <- meta_tbl[match(labels(dist_obj), meta_tbl$sample), , drop = FALSE]
  df$host <- factor(df$host)
  vegan::adonis2(dist_obj ~ host, data = df, by = "margin")
}

run_adonis_margin(dist_total, meta_for_plot)
run_adonis_margin(dist_endo, meta_for_plot)
run_adonis_margin(dist_nonendo, meta_for_plot)


# ----------------------------
# MASTER EXCEL SUMMARY WORKBOOK
# ----------------------------

if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}

excel_safe <- function(x) {
  x <- as.data.frame(x)
  for (nm in names(x)) {
    if (is.list(x[[nm]])) {
      x[[nm]] <- vapply(
        x[[nm]],
        function(z) paste(capture.output(str(z)), collapse = "; "),
        character(1)
      )
    }
  }
  x
}

read_optional_tsv <- function(path) {
  if (file.exists(path)) {
    readr::read_tsv(path, show_col_types = FALSE)
  } else {
    tibble::tibble()
  }
}

summarise_adonis <- function(dist_obj, meta_tbl, label) {
  df_perm <- meta_tbl[match(labels(dist_obj), meta_tbl$sample), , drop = FALSE]
  df_perm$host <- factor(df_perm$host)
  perm <- vegan::adonis2(dist_obj ~ host, data = df_perm)
  perm_df <- as.data.frame(perm)
  tibble::tibble(
    dataset = label,
    test = "PERMANOVA/adonis2",
    df = perm_df$Df[1],
    sum_of_squares = perm_df$SumOfSqs[1],
    r_squared = perm_df$R2[1],
    f_statistic = perm_df$F[1],
    p_value = perm_df$`Pr(>F)`[1]
  )
}

summarise_anosim <- function(dist_obj, meta_tbl, label) {
  df_ano <- meta_tbl[match(labels(dist_obj), meta_tbl$sample), , drop = FALSE]
  df_ano$host <- factor(df_ano$host)
  ano <- vegan::anosim(dist_obj, df_ano$host)
  tibble::tibble(
    dataset = label,
    test = "ANOSIM",
    statistic_R = as.numeric(ano$statistic),
    p_value = as.numeric(ano$signif),
    permutations = 999
  )
}

sample_overview <- meta_for_plot |>
  dplyr::left_join(valid_files |> dplyr::select(sample, report_path = path), by = "sample") |>
  dplyr::arrange(host, sample)

abundance_summary <- tibble::tibble(
  dataset = c("total", "endo", "nonendo"),
  samples = c(nrow(total_mat), nrow(endo_mat), nrow(nonendo_mat)),
  taxa = c(ncol(total_mat), ncol(endo_mat), ncol(nonendo_mat)),
  total_reads = c(sum(total_mat), sum(endo_mat), sum(nonendo_mat)),
  median_reads_per_sample = c(
    median(rowSums(total_mat)),
    median(rowSums(endo_mat)),
    median(rowSums(nonendo_mat))
  )
)

alpha_by_sample <- dplyr::bind_rows(
  alpha_total |> dplyr::mutate(dataset = "total"),
  alpha_endo |> dplyr::mutate(dataset = "endo"),
  alpha_nonendo |> dplyr::mutate(dataset = "nonendo")
) |>
  dplyr::left_join(meta_for_plot, by = "sample") |>
  dplyr::select(dataset, sample, host, richness, shannon, simpson)

alpha_by_host <- alpha_by_sample |>
  dplyr::group_by(dataset, host) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_richness = mean(richness, na.rm = TRUE),
    median_richness = median(richness, na.rm = TRUE),
    mean_shannon = mean(shannon, na.rm = TRUE),
    median_shannon = median(shannon, na.rm = TRUE),
    mean_simpson = mean(simpson, na.rm = TRUE),
    median_simpson = median(simpson, na.rm = TRUE),
    .groups = "drop"
  )

pcoa_variance <- dplyr::bind_rows(
  read_optional_tsv(file.path(pcoa_dirs$data, "pcoa_total_variance.tsv")) |> dplyr::mutate(dataset = "total"),
  read_optional_tsv(file.path(pcoa_dirs$data, "pcoa_endo_variance.tsv")) |> dplyr::mutate(dataset = "endo"),
  read_optional_tsv(file.path(pcoa_dirs$data, "pcoa_nonendo_variance.tsv")) |> dplyr::mutate(dataset = "nonendo")
) |>
  dplyr::select(dataset, dplyr::everything())

pcoa_envfit <- dplyr::bind_rows(
  read_optional_tsv(file.path(pcoa_dirs$results, "pcoa_total_envfit.tsv")) |> dplyr::mutate(dataset = "total"),
  read_optional_tsv(file.path(pcoa_dirs$results, "pcoa_endo_envfit.tsv")) |> dplyr::mutate(dataset = "endo"),
  read_optional_tsv(file.path(pcoa_dirs$results, "pcoa_nonendo_envfit.tsv")) |> dplyr::mutate(dataset = "nonendo")
) |>
  dplyr::select(dataset, dplyr::everything())

permanova_summary <- dplyr::bind_rows(
  summarise_adonis(dist_total, meta_for_plot, "total"),
  summarise_adonis(dist_endo, meta_for_plot, "endo"),
  summarise_adonis(dist_nonendo, meta_for_plot, "nonendo")
)

anosim_summary <- dplyr::bind_rows(
  summarise_anosim(dist_total, meta_for_plot, "total"),
  summarise_anosim(dist_endo, meta_for_plot, "endo"),
  summarise_anosim(dist_nonendo, meta_for_plot, "nonendo")
)

deseq_lrt_summary <- res_lrt_df |>
  dplyr::arrange(padj) |>
  dplyr::mutate(rank_by_padj = dplyr::row_number()) |>
  dplyr::select(rank_by_padj, taxon, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

coverage_summary_for_excel <- if (exists("coverage_stats_summary")) {
  coverage_stats_summary
} else {
  tibble::tibble()
}

key_outputs <- tibble::tibble(
  analysis = c(
    "abundance_matrices",
    "alpha_diversity",
    "pcoa_ordination",
    "beta_diversity_tests",
    "permdisp",
    "heatmaps",
    "differential_abundance",
    "family_endogenization",
    "coverage_analysis"
  ),
  data_folder = c(
    abundance_dirs$data,
    alpha_dirs$data,
    pcoa_dirs$data,
    beta_dirs$data,
    permdisp_dirs$data,
    heatmap_dirs$data,
    deseq_dirs$data,
    family_dirs$data,
    coverage_dirs$data
  ),
  results_folder = c(
    abundance_dirs$results,
    alpha_dirs$results,
    pcoa_dirs$results,
    beta_dirs$results,
    permdisp_dirs$results,
    heatmap_dirs$results,
    deseq_dirs$results,
    family_dirs$results,
    coverage_dirs$results
  ),
  figures_folder = c(
    abundance_dirs$figures,
    alpha_dirs$figures,
    pcoa_dirs$figures,
    beta_dirs$figures,
    permdisp_dirs$figures,
    heatmap_dirs$figures,
    deseq_dirs$figures,
    family_dirs$figures,
    coverage_dirs$figures
  )
)

master_readme <- tibble::tibble(
  field = c(
    "generated_on",
    "input_metadata_file",
    "input_summary_file",
    "kraken_report_directories",
    "master_workbook_description"
  ),
  value = c(
    as.character(Sys.time()),
    metadata_file,
    summary_file,
    paste(c(alb_dir, cpip_dir, cmel_dir), collapse = "; "),
    "Summary workbook for abundance matrices, alpha diversity, ordination, beta-diversity tests, PERMDISP, differential abundance, family endogenization, and coverage analyses."
  )
)

master_sheets <- list(
  README = excel_safe(master_readme),
  folder_index = excel_safe(key_outputs),
  sample_overview = excel_safe(sample_overview),
  abundance_summary = excel_safe(abundance_summary),
  alpha_by_sample = excel_safe(alpha_by_sample),
  alpha_by_host = excel_safe(alpha_by_host),
  pcoa_variance = excel_safe(pcoa_variance),
  pcoa_envfit = excel_safe(pcoa_envfit),
  permanova = excel_safe(permanova_summary),
  anosim = excel_safe(anosim_summary),
  permdisp = excel_safe(summary_rows),
  deseq_host_LRT = excel_safe(deseq_lrt_summary),
  family_totals = excel_safe(family_totals),
  family_fisher = excel_safe(stats_list),
  family_collapsed = excel_safe(collapsed),
  coverage_wilcox = excel_safe(coverage_summary_for_excel)
)

master_excel_path <- file.path(outdir, "master_analysis_summary.xlsx")
writexl::write_xlsx(master_sheets, path = master_excel_path)

cat("\nMaster Excel summary written to:", normalizePath(master_excel_path), "\n")

cat("\nAll analyses complete. Organized outputs written under:", normalizePath(outdir), "\n")


# ============================================================
# Within-host spatiotemporal PERMANOVA
# ============================================================

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------

collection_metadata_file <- "sample_collection_metadata.csv"
outdir <- "diversity_results/11_within_host_spatiotemporal"

permutations_n <- 499

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "data"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "logs"), showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

required_packages <- c("tidyverse", "vegan")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

# ------------------------------------------------------------
# Required objects from the main diversity script
# ------------------------------------------------------------

required_objects <- c("dist_total", "dist_endo", "dist_nonendo")

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), envir = .GlobalEnv)]

if (length(missing_objects) > 0) {
  stop(
    "This add-on script must be run after the main diversity script in the same R session.\n",
    "Missing objects: ", paste(missing_objects, collapse = ", "), "\n\n",
    "Run the main script first, then run this script."
  )
}

# ------------------------------------------------------------
# Load and standardize collection metadata
# ------------------------------------------------------------

collection_meta_raw <- readr::read_csv(collection_metadata_file, show_col_types = FALSE)

names(collection_meta_raw) <- names(collection_meta_raw) |>
  stringr::str_trim() |>
  stringr::str_replace_all("\\s+", "_") |>
  stringr::str_replace_all("\\.", "_") |>
  tolower()

find_col <- function(df, candidates, required = TRUE, label = "column") {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) {
    if (required) {
      stop(
        "Could not find ", label, ". Tried: ",
        paste(candidates, collapse = ", "),
        "\nAvailable columns are: ",
        paste(names(df), collapse = ", ")
      )
    } else {
      return(NA_character_)
    }
  }
  hit[1]
}

sample_col <- find_col(collection_meta_raw, c("sample", "sample_number", "sample_id"), TRUE, "sample column")
host_col <- find_col(collection_meta_raw, c("host", "species"), TRUE, "host column")
site_col <- find_col(collection_meta_raw, c("site", "collection_site"), FALSE, "site column")
county_col <- find_col(collection_meta_raw, c("county"), FALSE, "county column")
town_col <- find_col(collection_meta_raw, c("town", "city"), FALSE, "town column")
date_col <- find_col(collection_meta_raw, c("date", "collection_date"), FALSE, "date column")
lat_col <- find_col(collection_meta_raw, c("latitude", "lat"), FALSE, "latitude column")
lon_col <- find_col(collection_meta_raw, c("longitude", "long", "lon"), FALSE, "longitude column")

collection_meta <- collection_meta_raw |>
  dplyr::transmute(
    sample = as.character(.data[[sample_col]]),
    host = as.character(.data[[host_col]]),
    site = if (!is.na(site_col)) as.character(.data[[site_col]]) else NA_character_,
    county = if (!is.na(county_col)) as.character(.data[[county_col]]) else NA_character_,
    town = if (!is.na(town_col)) as.character(.data[[town_col]]) else NA_character_,
    collection_date_raw = if (!is.na(date_col)) as.character(.data[[date_col]]) else NA_character_,
    latitude = if (!is.na(lat_col)) suppressWarnings(as.numeric(.data[[lat_col]])) else NA_real_,
    longitude = if (!is.na(lon_col)) suppressWarnings(as.numeric(.data[[lon_col]])) else NA_real_
  ) |>
  dplyr::mutate(
    collection_date = dplyr::coalesce(
      suppressWarnings(as.Date(collection_date_raw, format = "%Y-%m-%d")),
      suppressWarnings(as.Date(collection_date_raw, format = "%m/%d/%Y")),
      suppressWarnings(as.Date(collection_date_raw, format = "%d-%b-%y")),
      suppressWarnings(as.Date(collection_date_raw, format = "%d-%b-%Y"))
    ),
    collection_month_num = suppressWarnings(as.integer(format(collection_date, "%m"))),
    collection_month = ifelse(
      is.na(collection_month_num),
      NA_character_,
      month.name[collection_month_num]
    ),
    collection_month = factor(collection_month, levels = month.name, ordered = TRUE),
    host = factor(
      host,
      levels = c("Aedes albopictus", "Culex pipiens", "Culiseta melanura")
    ),
    site = factor(site),
    county = factor(county),
    town = factor(town)
  ) |>
  dplyr::distinct(sample, .keep_all = TRUE)

readr::write_tsv(
  collection_meta,
  file.path(outdir, "data", "collection_metadata_standardized.tsv")
)

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

tidy_adonis2 <- function(x) {
  out <- as.data.frame(x)
  out$term <- rownames(out)
  rownames(out) <- NULL
  
  names(out) <- names(out) |>
    stringr::str_replace_all("SumOfSqs", "sum_of_squares") |>
    stringr::str_replace_all("R2", "r_squared") |>
    stringr::str_replace_all("^F$", "f_statistic") |>
    stringr::str_replace_all("Pr\\(>F\\)", "p_value")
  
  out |>
    dplyr::select(term, dplyr::everything())
}

get_dist_samples <- function(dist_obj) {
  labels(dist_obj)
}

subset_dist <- function(dist_obj, samples) {
  dmat <- as.matrix(dist_obj)
  samples <- intersect(samples, rownames(dmat))
  if (length(samples) < 3) return(NULL)
  stats::as.dist(dmat[samples, samples, drop = FALSE])
}

has_enough_factor_levels <- function(df, variable) {
  if (!variable %in% names(df)) return(FALSE)
  x <- droplevels(as.factor(df[[variable]]))
  x <- x[!is.na(x)]
  if (length(x) < 3) return(FALSE)
  nlevels(x) >= 2
}

has_enough_numeric_values <- function(df, variable) {
  if (!variable %in% names(df)) return(FALSE)
  x <- suppressWarnings(as.numeric(df[[variable]]))
  sum(!is.na(x)) >= 4 && length(unique(x[!is.na(x)])) >= 2
}

model_is_reasonable <- function(df, model_name) {
  n <- nrow(df)
  
  if (n < 5) {
    return(list(ok = FALSE, reason = paste0("too few samples: n = ", n)))
  }
  
  if (model_name == "county") {
    if (!has_enough_factor_levels(df, "county")) {
      return(list(ok = FALSE, reason = "county has fewer than two observed levels"))
    }
    if (nlevels(droplevels(df$county)) >= n - 1) {
      return(list(ok = FALSE, reason = "county has too many levels for available samples"))
    }
  }
  
  if (model_name == "site") {
    if (!has_enough_factor_levels(df, "site")) {
      return(list(ok = FALSE, reason = "site has fewer than two observed levels"))
    }
    if (nlevels(droplevels(df$site)) >= n - 1) {
      return(list(ok = FALSE, reason = "site has too many levels for available samples"))
    }
  }
  
  if (model_name == "collection_month") {
    if (!has_enough_factor_levels(df, "collection_month")) {
      return(list(ok = FALSE, reason = "collection_month has fewer than two observed levels"))
    }
    if (nlevels(droplevels(df$collection_month)) >= n - 1) {
      return(list(ok = FALSE, reason = "collection_month has too many levels for available samples"))
    }
  }
  
  if (model_name == "county_plus_month") {
    if (!has_enough_factor_levels(df, "county")) {
      return(list(ok = FALSE, reason = "county has fewer than two observed levels"))
    }
    if (!has_enough_factor_levels(df, "collection_month")) {
      return(list(ok = FALSE, reason = "collection_month has fewer than two observed levels"))
    }
    
    df_model <- df |>
      dplyr::filter(!is.na(county), !is.na(collection_month)) |>
      dplyr::mutate(
        county = droplevels(county),
        collection_month = droplevels(collection_month)
      )
    
    approx_parameters <- nlevels(df_model$county) + nlevels(df_model$collection_month)
    if (nrow(df_model) <= approx_parameters + 2) {
      return(list(ok = FALSE, reason = "county + month model is over-parameterized for available samples"))
    }
  }
  
  if (model_name == "lat_long_month_numeric") {
    needed <- c("latitude", "longitude", "collection_month_num")
    if (!all(vapply(needed, function(v) has_enough_numeric_values(df, v), logical(1)))) {
      return(list(ok = FALSE, reason = "latitude, longitude, or numeric month has insufficient numeric variation"))
    }
    
    df_model <- df |>
      dplyr::filter(!is.na(latitude), !is.na(longitude), !is.na(collection_month_num))
    
    if (nrow(df_model) < 8) {
      return(list(ok = FALSE, reason = "too few complete samples for latitude + longitude + month model"))
    }
  }
  
  list(ok = TRUE, reason = "ok")
}

run_safe_adonis <- function(dist_obj, meta_tbl, dataset, host_name, model_name, formula_text) {
  host_safe <- gsub("[^A-Za-z0-9]+", "_", host_name)
  result_file <- file.path(
    outdir,
    "results",
    paste0("permanova_", dataset, "_", host_safe, "_", model_name, ".txt")
  )
  
  samples_in_dist <- get_dist_samples(dist_obj)
  
  df <- meta_tbl |>
    dplyr::filter(sample %in% samples_in_dist, as.character(host) == host_name) |>
    dplyr::arrange(match(sample, samples_in_dist))
  
  reason <- model_is_reasonable(df, model_name)
  
  if (!reason$ok) {
    writeLines(
      c(
        paste0("Dataset: ", dataset),
        paste0("Host: ", host_name),
        paste0("Model: ", model_name),
        paste0("Formula: ", formula_text),
        paste0("Status: skipped"),
        paste0("Reason: ", reason$reason)
      ),
      con = result_file
    )
    
    return(tibble::tibble(
      dataset = dataset,
      host = host_name,
      model = model_name,
      term = NA_character_,
      df = NA_real_,
      sum_of_squares = NA_real_,
      r_squared = NA_real_,
      f_statistic = NA_real_,
      p_value = NA_real_,
      status = "skipped",
      reason = reason$reason,
      n_samples = nrow(df)
    ))
  }
  
  formula_vars <- all.vars(stats::as.formula(formula_text))[-1]
  
  df_model <- df |>
    dplyr::filter(dplyr::if_all(dplyr::all_of(formula_vars), ~ !is.na(.x))) |>
    dplyr::mutate(
      county = droplevels(county),
      site = droplevels(site),
      collection_month = droplevels(collection_month)
    )
  
  dist_sub <- subset_dist(dist_obj, df_model$sample)
  
  if (is.null(dist_sub) || length(labels(dist_sub)) < 5) {
    msg <- "too few samples after filtering complete metadata"
    
    writeLines(
      c(
        paste0("Dataset: ", dataset),
        paste0("Host: ", host_name),
        paste0("Model: ", model_name),
        paste0("Formula: ", formula_text),
        paste0("Status: skipped"),
        paste0("Reason: ", msg)
      ),
      con = result_file
    )
    
    return(tibble::tibble(
      dataset = dataset,
      host = host_name,
      model = model_name,
      term = NA_character_,
      df = NA_real_,
      sum_of_squares = NA_real_,
      r_squared = NA_real_,
      f_statistic = NA_real_,
      p_value = NA_real_,
      status = "skipped",
      reason = msg,
      n_samples = nrow(df_model)
    ))
  }
  
  df_model <- df_model[match(labels(dist_sub), df_model$sample), , drop = FALSE]
  
  model_result <- tryCatch(
    {
      vegan::adonis2(
        stats::as.formula(formula_text),
        data = df_model,
        permutations = permutations_n,
        by = "margin"
      )
    },
    error = function(e) e
  )
  
  if (inherits(model_result, "error")) {
    msg <- conditionMessage(model_result)
    
    writeLines(
      c(
        paste0("Dataset: ", dataset),
        paste0("Host: ", host_name),
        paste0("Model: ", model_name),
        paste0("Formula: ", formula_text),
        paste0("Status: error"),
        paste0("Reason: ", msg)
      ),
      con = result_file
    )
    
    return(tibble::tibble(
      dataset = dataset,
      host = host_name,
      model = model_name,
      term = NA_character_,
      df = NA_real_,
      sum_of_squares = NA_real_,
      r_squared = NA_real_,
      f_statistic = NA_real_,
      p_value = NA_real_,
      status = "error",
      reason = msg,
      n_samples = nrow(df_model)
    ))
  }
  
  capture.output(model_result, file = result_file)
  
  tidy_result <- tidy_adonis2(model_result)
  
  tidy_result |>
    dplyr::mutate(
      dataset = dataset,
      host = host_name,
      model = model_name,
      status = "ok",
      reason = NA_character_,
      n_samples = nrow(df_model),
      .before = 1
    ) |>
    dplyr::rename(df = Df) |>
    dplyr::select(
      dataset,
      host,
      model,
      term,
      df,
      dplyr::any_of(c("sum_of_squares", "r_squared", "f_statistic", "p_value")),
      status,
      reason,
      n_samples
    )
}

# ------------------------------------------------------------
# Run within-host models
# ------------------------------------------------------------

dist_list <- list(
  total = dist_total,
  endo = dist_endo,
  nonendo = dist_nonendo
)

hosts_to_test <- collection_meta |>
  dplyr::filter(!is.na(host)) |>
  dplyr::pull(host) |>
  as.character() |>
  unique()

model_list <- tibble::tribble(
  ~model_name, ~formula_text,
  "county", "dist_sub ~ county",
  "site", "dist_sub ~ site",
  "collection_month", "dist_sub ~ collection_month",
  "county_plus_month", "dist_sub ~ county + collection_month",
  "lat_long_month_numeric", "dist_sub ~ latitude + longitude + collection_month_num"
)

all_results <- list()

counter <- 1

for (dataset_name in names(dist_list)) {
  message("Dataset: ", dataset_name)
  
  for (host_name in hosts_to_test) {
    message("  Host: ", host_name)
    
    for (i in seq_len(nrow(model_list))) {
      model_name <- model_list$model_name[i]
      formula_text <- model_list$formula_text[i]
      
      message("    Model: ", model_name)
      
      all_results[[counter]] <- run_safe_adonis(
        dist_obj = dist_list[[dataset_name]],
        meta_tbl = collection_meta,
        dataset = dataset_name,
        host_name = host_name,
        model_name = model_name,
        formula_text = formula_text
      )
      
      counter <- counter + 1
      gc(verbose = FALSE)
    }
  }
}

summary_tbl <- dplyr::bind_rows(all_results)

readr::write_tsv(
  summary_tbl,
  file.path(outdir, "results", "within_host_spatiotemporal_permanova_summary.tsv")
)

significant_tbl <- summary_tbl |>
  dplyr::filter(
    status == "ok",
    !is.na(p_value),
    p_value < 0.05,
    !term %in% c("Residual", "Total")
  )

readr::write_tsv(
  significant_tbl,
  file.path(outdir, "results", "within_host_spatiotemporal_permanova_significant_only.tsv")
)

# ------------------------------------------------------------
# Save run information
# ------------------------------------------------------------

sink(file.path(outdir, "logs", "within_host_spatiotemporal_session_info.txt"))
cat("Within-host spatiotemporal PERMANOVA completed:\n")
cat(as.character(Sys.time()), "\n\n")
cat("Permutations:", permutations_n, "\n\n")
cat("Input metadata:", collection_metadata_file, "\n\n")
cat("Hosts tested:\n")
print(hosts_to_test)
cat("\nModels tested:\n")
print(model_list)
cat("\nSummary of model status:\n")
print(table(summary_tbl$status, useNA = "ifany"))
cat("\nSession info:\n")
sessionInfo()
sink()

cat("\nWithin-host spatiotemporal PERMANOVA complete.\n")
cat("Main summary table:\n")
cat(file.path(outdir, "results", "within_host_spatiotemporal_permanova_summary.tsv"), "\n")
cat("\nSignificant-only table:\n")
cat(file.path(outdir, "results", "within_host_spatiotemporal_permanova_significant_only.tsv"), "\n")


library(tidyverse)
library(lubridate)

# Choose which dataset you want to test:
# Options: "total", "endo", "nonendo"
dataset <- "total"

pcoa_file <- file.path(
  "diversity_results",
  "03_pcoa_ordination",
  "data",
  paste0("pcoa_", dataset, "_host.tsv")
)

metadata_file <- "sample_collection_metadata.csv"

pcoa_scores <- readr::read_tsv(pcoa_file, show_col_types = FALSE)

collection_meta <- readr::read_csv(metadata_file, show_col_types = FALSE) |>
  janitor::clean_names() |>
  dplyr::rename(
    sample = sample,
    host_metadata = host
  ) |>
  dplyr::mutate(
    date = lubridate::mdy(date),
    collection_month = lubridate::month(date, label = TRUE, abbr = FALSE),
    collection_month_num = lubridate::month(date),
    county = factor(county),
    site = factor(site)
  )

pcoa_df <- pcoa_scores |>
  dplyr::left_join(collection_meta, by = "sample") |>
  dplyr::mutate(
    host = factor(host),
    county = factor(county),
    site = factor(site),
    collection_month = factor(collection_month)
  )

# Check that the join worked
dplyr::glimpse(pcoa_df)

lm_pcoa1 <- lm(
  PCoA1 ~ host + county + collection_month + latitude + longitude,
  data = pcoa_df
)

summary(lm_pcoa1)
anova(lm_pcoa1)


lm_pcoa1 <- lm(PCoA1 ~ host + county + collection_month + latitude + longitude, data = pcoa_df)
lm_pcoa2 <- lm(PCoA2 ~ host + county + collection_month + latitude + longitude, data = pcoa_df)
lm_pcoa3 <- lm(PCoA3 ~ host + county + collection_month + latitude + longitude, data = pcoa_df)

summary(lm_pcoa1)
summary(lm_pcoa2)
summary(lm_pcoa3)

anova(lm_pcoa1)
anova(lm_pcoa2)
anova(lm_pcoa3)


# ------------------------------------------------------------
# Save PCoA axis model results to Excel
# ------------------------------------------------------------

#required_packages <- c("broom", "openxlsx", "dplyr", "purrr", "tibble", "readr")

#for (pkg in required_packages) {
#  if (!requireNamespace(pkg, quietly = TRUE)) {
#    install.packages(pkg)
#  }
#}

library(broom)
library(openxlsx)
library(dplyr)
library(purrr)
library(tibble)
library(readr)

# Put models into a named list
pcoa_models <- list(
  PCoA1 = lm_pcoa1,
  PCoA2 = lm_pcoa2,
  PCoA3 = lm_pcoa3
)

# ----------------------------
# Coefficient-level results
# ----------------------------

coef_results <- purrr::map_dfr(
  names(pcoa_models),
  function(axis_name) {
    broom::tidy(pcoa_models[[axis_name]]) |>
      dplyr::mutate(
        axis = axis_name,
        result_type = "coefficient"
      ) |>
      dplyr::select(
        axis,
        result_type,
        term,
        estimate,
        std.error,
        statistic,
        p.value
      )
  }
)

# ----------------------------
# ANOVA term-level results
# ----------------------------

anova_results <- purrr::map_dfr(
  names(pcoa_models),
  function(axis_name) {
    anova_tbl <- as.data.frame(anova(pcoa_models[[axis_name]]))
    
    anova_tbl |>
      tibble::rownames_to_column("term") |>
      dplyr::mutate(
        axis = axis_name,
        result_type = "anova"
      ) |>
      dplyr::rename(
        df = Df,
        sum_sq = `Sum Sq`,
        mean_sq = `Mean Sq`,
        f_value = `F value`,
        p_value = `Pr(>F)`
      ) |>
      dplyr::select(
        axis,
        result_type,
        term,
        df,
        sum_sq,
        mean_sq,
        f_value,
        p_value
      )
  }
)

# ----------------------------
# Model-level summaries
# ----------------------------

model_summary <- purrr::map_dfr(
  names(pcoa_models),
  function(axis_name) {
    glance_tbl <- broom::glance(pcoa_models[[axis_name]])
    
    tibble::tibble(
      axis = axis_name,
      r_squared = glance_tbl$r.squared,
      adjusted_r_squared = glance_tbl$adj.r.squared,
      sigma = glance_tbl$sigma,
      model_f_statistic = glance_tbl$statistic,
      model_p_value = glance_tbl$p.value,
      df_model = glance_tbl$df,
      df_residual = glance_tbl$df.residual,
      AIC = glance_tbl$AIC,
      BIC = glance_tbl$BIC,
      deviance = glance_tbl$deviance
    )
  }
)

# ----------------------------
# Combined main results sheet
# ----------------------------

coef_results_for_main <- coef_results |>
  dplyr::rename(
    std_error = std.error,
    t_value = statistic,
    p_value = p.value
  ) |>
  dplyr::mutate(
    df = NA_real_,
    sum_sq = NA_real_,
    mean_sq = NA_real_,
    f_value = NA_real_
  ) |>
  dplyr::select(
    axis,
    result_type,
    term,
    estimate,
    std_error,
    t_value,
    df,
    sum_sq,
    mean_sq,
    f_value,
    p_value
  )

anova_results_for_main <- anova_results |>
  dplyr::mutate(
    estimate = NA_real_,
    std_error = NA_real_,
    t_value = NA_real_
  ) |>
  dplyr::select(
    axis,
    result_type,
    term,
    estimate,
    std_error,
    t_value,
    df,
    sum_sq,
    mean_sq,
    f_value,
    p_value
  )

pcoa_axes_results <- dplyr::bind_rows(
  coef_results_for_main,
  anova_results_for_main
) |>
  dplyr::mutate(
    significance = dplyr::case_when(
      is.na(p_value) ~ "",
      p_value <= 0.001 ~ "***",
      p_value <= 0.01 ~ "**",
      p_value <= 0.05 ~ "*",
      p_value <= 0.1 ~ ".",
      TRUE ~ "ns"
    )
  )

# ----------------------------
# Optional: save CSV copies too
# ----------------------------

readr::write_csv(pcoa_axes_results, "pcoa_axes_results.csv")
readr::write_csv(coef_results, "pcoa_axes_coefficients.csv")
readr::write_csv(anova_results, "pcoa_axes_anova.csv")
readr::write_csv(model_summary, "pcoa_axes_model_summary.csv")

# ----------------------------
# Create Excel workbook
# ----------------------------

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "pcoa_axes_results")
openxlsx::addWorksheet(wb, "coefficients")
openxlsx::addWorksheet(wb, "anova")
openxlsx::addWorksheet(wb, "model_summary")

openxlsx::writeData(wb, "pcoa_axes_results", pcoa_axes_results)
openxlsx::writeData(wb, "coefficients", coef_results)
openxlsx::writeData(wb, "anova", anova_results)
openxlsx::writeData(wb, "model_summary", model_summary)

# Add basic formatting
header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom"
)

for (sheet in names(wb)) {
  openxlsx::addStyle(
    wb,
    sheet = sheet,
    style = header_style,
    rows = 1,
    cols = 1:ncol(openxlsx::readWorkbook(wb, sheet = sheet)),
    gridExpand = TRUE
  )
  
  openxlsx::freezePane(wb, sheet = sheet, firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet = sheet, cols = 1:50, widths = "auto")
}

openxlsx::saveWorkbook(
  wb,
  file = "pcoa_axes_results.xlsx",
  overwrite = TRUE
)

cat("\nSaved PCoA axis model results to:\n")
cat("pcoa_axes_results.xlsx\n")
cat("pcoa_axes_results.csv\n")
cat("pcoa_axes_coefficients.csv\n")
cat("pcoa_axes_anova.csv\n")
cat("pcoa_axes_model_summary.csv\n")