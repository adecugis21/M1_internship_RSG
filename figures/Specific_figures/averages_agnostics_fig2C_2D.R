#!/usr/bin/env Rscript
# ---------------------------------------------------------
# Average profiles + correlations for H3K4ME3
# ---------------------------------------------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))

suppressPackageStartupMessages({
  library(tidyCoverage)
  library(tidySummarizedExperiment)
  library(rtracklayer)
  library(plyranges)
  library(purrr)
  library(ggplot2)
  library(BiocParallel)
  library(dplyr)
  library(tidyr)
})

# ---------------------------------------------------------
# SETTINGS
# ---------------------------------------------------------
kb_around <- 5000   # distance around TSSs
species   <- "cerevisiae"

bw_folder <- file.path(
  "/pasteur/helix/scratch/adecugis/all_momics/momics/output/predictions/figures_methylation_multispecies_no_maxpool",
  species
)
bw_folder_single <- file.path(
  "/pasteur/helix/scratch/adecugis/all_momics/momics/output/predictions/figures_one_track_multispecies_no_maxpool",
  species
)


# ---------------------------------------------------------
# FILES
# ---------------------------------------------------------
bw_files <- c(
  ground_truth_H3K4ME3 = file.path(
    "/pasteur/helix/scratch/adecugis/all_momics/momics/output/a_extracting_bw", species,
    paste0(species, "_H3K4ME3_rescaled.bw")
  ),
  prediction_H3K4ME3_single_species = file.path(
    bw_folder_single,
    paste0(species, "_pred_XV_H3K4ME3.bw")
  ),
  prediction_H3K4ME3_agnostic = file.path(
    bw_folder,
    paste0(species, "_pred_XV_H3K4ME3.bw")
  )
)


tracks <- BigWigFileList(bw_files)

# ---------------------------------------------------------
# FEATURES — pick TSS BED depending on species
# ---------------------------------------------------------
tss_bed_files <- list(
  cerevisiae = system.file("extdata", "TSSs.bed", package = "tidyCoverage"),
  pombe      = "/pasteur/appa/scratch/adecugis/all_momics/momics/scripts/Average_Profiles/input_out/pombe_TSS.bed",
  glabrata   = "/pasteur/appa/scratch/adecugis/all_momics/momics/scripts/Average_Profiles/input_out/Candida_glabrata_TSSs.bed"
)
features <- list(TSSs = import(tss_bed_files[[species]]))

# ---------------------------------------------------------
# CoverageExperiment
# ---------------------------------------------------------
CE <- CoverageExperiment(
  tracks,
  features,
  width        = kb_around,
  ignore.strand = FALSE,
  BPPARAM      = BiocParallel::SerialParam()
)
CE  <- dplyr::filter(CE, features == "TSSs")
agg <- aggregate(CE)
print(head(as_tibble(agg)))
stopifnot(nrow(agg) > 0)

# ---------------------------------------------------------
# Wide format for correlations
# FIX: explicit dplyr:: / tidyr:: to avoid masking conflicts
# ---------------------------------------------------------
df_wide <- agg %>%
  as_tibble() %>%
  dplyr::select(coord, track, mean) %>%
  tidyr::pivot_wider(names_from = track, values_from = mean)

# ---------------------------------------------------------
# Track labels and color setup (define before use)
# ---------------------------------------------------------
track_labels <- c(
  ground_truth_H3K4ME3              = "H3K4ME3 Ground Truth",
  prediction_H3K4ME3_single_species = "H3K4ME3 \u2013 Single-species model",
  prediction_H3K4ME3_agnostic       = "H3K4ME3 \u2013 Cross-species model"
)

cols <- RColorBrewer::brewer.pal(9, "Set1")
names(cols) <- track_labels

role_colors <- c(
  "Ground Truth"   = "#1B9E77",
  "Single-species" = "#D95F02",
  "Cross-species"  = "#7570B3"
)

# Convert df to long format for plotting
df <- df_wide %>%
  tidyr::pivot_longer(
    cols = -coord,
    names_to = "track",
    values_to = "coverage"
  ) %>%
  dplyr::mutate(
    position = coord,
    track = factor(track, levels = names(track_labels), labels = track_labels),
    track_role = dplyr::case_when(
      grepl("Ground Truth", as.character(track), fixed = TRUE) ~ "Ground Truth",
      grepl("Single-species", as.character(track), fixed = TRUE) ~ "Single-species",
      grepl("Cross-species|Multispecies", as.character(track), ignore.case = TRUE) ~ "Cross-species",
      TRUE ~ NA_character_
    ),
    mark = "H3K4ME3"
  )

# ---------------------------------------------------------
# Compute correlations
# ---------------------------------------------------------
corr_tbl <- tibble(
  mark = "H3K4ME3",
  pearson_single = cor(
    df_wide$ground_truth_H3K4ME3,
    df_wide$prediction_H3K4ME3_single_species,
    use = "complete.obs"
  ),
  spearman_single = cor(
    df_wide$ground_truth_H3K4ME3,
    df_wide$prediction_H3K4ME3_single_species,
    method = "spearman",
    use = "complete.obs"
  ),
  pearson_agnostic = cor(
    df_wide$ground_truth_H3K4ME3,
    df_wide$prediction_H3K4ME3_agnostic,
    use = "complete.obs"
  ),
  spearman_agnostic = cor(
    df_wide$ground_truth_H3K4ME3,
    df_wide$prediction_H3K4ME3_agnostic,
    method = "spearman",
    use = "complete.obs"
  )
)

# ---------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------
make_panel <- function(df_sub, ylim, title = NULL, subtitle = NULL, corr_text = NULL) {
  x_range <- range(df_sub$position, na.rm = TRUE)
  x_breaks <- unique(c(x_range[1], 0, x_range[2]))
  x_labels <- ifelse(x_breaks == 0, "TSS", sprintf("%.1f kb", x_breaks / 1000))

  p <- ggplot(df_sub, aes(x = position, y = coverage, colour = track_role)) +
    stat_summary(
      fun       = mean,
      geom      = "line",
      linewidth = 0.7,
      alpha     = 0.85
    ) +
    scale_colour_manual(values = role_colors) +
    scale_x_continuous(
      limits = x_range,
      breaks = x_breaks,
      labels = x_labels,
      expand = c(0.02, 0)
    ) +
    coord_cartesian(ylim = ylim) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = "Distance from TSS (bp)",
      y        = "Normalised signal coverage around TSS of chromosome XV forward strand",
      colour   = "Track"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle   = element_text(hjust = 0.5, size = 12),
      legend.position = "none",
      axis.text       = element_text(colour = "black")
    )

  p
}

make_mark_section <- function() {
  df_mark <- dplyr::filter(df, mark == "H3K4ME3")

  gt_label <- track_labels[["ground_truth_H3K4ME3"]]
  single_label <- track_labels[["prediction_H3K4ME3_single_species"]]
  multi_label <- track_labels[["prediction_H3K4ME3_agnostic"]]

  df_overlay <- dplyr::filter(df_mark, track %in% c(gt_label, single_label, multi_label))
  df_pred_only <- dplyr::filter(df_mark, track %in% c(single_label, multi_label))

  overlay_ylim <- range(df_overlay$coverage, na.rm = TRUE)
  pred_ylim    <- range(df_pred_only$coverage, na.rm = TRUE)

  corr_single <- sprintf(
    "Pearson r = %.3f  |  Spearman ρ = %.3f",
    corr_tbl$pearson_single,
    corr_tbl$spearman_single
  )
  corr_multi <- sprintf(
    "Pearson r = %.3f  |  Spearman ρ = %.3f",
    corr_tbl$pearson_agnostic,
    corr_tbl$spearman_agnostic
  )

  subtitle_text <- paste(
    "single-species correlations:", corr_single,
    "\ncross-species correlations:", corr_multi,
    "\n\nGround Truth (green) | Single-track (orange) | Multi-modal (blue)"
  )

  list(
    make_panel(
      df_overlay,
      ylim = overlay_ylim,
      title = paste("S. cerevisiae - H3K4ME3 signal around TSSs: ground truth vs predictions"),
      subtitle = subtitle_text
    ),
    make_panel(
      df_pred_only,
      ylim = pred_ylim,
      title = paste("Single-track  vs multi-modal predictions"),
    )
  )
}

plot_list <- make_mark_section()

# ---------------------------------------------------------
# Save first plot
# ---------------------------------------------------------
outdir <- "/pasteur/helix/scratch/adecugis/all_momics/momics/output/R/figures_rapport"
png_out <- file.path(
  outdir,
  paste0("coverage_", species, "_H3K4ME3_2panels_sharedscale.png")
)

png(filename = png_out, width = 8, height = 12, units = "in", res = 400)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 2, ncol = 1)))

for (i in seq_along(plot_list)) {
  print(plot_list[[i]], vp = grid::viewport(layout.pos.row = i, layout.pos.col = 1))
}

dev.off()
message("\u2705 Plot saved to: ", png_out)