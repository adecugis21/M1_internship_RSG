#!/usr/bin/env Rscript
options(repos = c(CRAN = "https://cloud.r-project.org"))

library(tidyCoverage)
library(tidySummarizedExperiment)
library(rtracklayer)
library(plyranges)
library(ggplot2)
library(BiocParallel)
library(dplyr)
library(magick)

# Add this helper function near the top of the script (after library calls)
get_track_colors <- function(n) {
  base_colors <- c(
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666",
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
    "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
    "#9C755F", "#BAB0AC"
  )
  if (n <= length(base_colors)) return(base_colors[1:n])
  return(colorRampPalette(base_colors)(n))
}

# =========================================================
# USER SETTINGS
# =========================================================
tag       <- "May_26_analysis"
kb_around <- 5000

PRED_DIRS <- list(
  concatenated = list(
    path             = "/pasteur/helix/scratch/adecugis/all_momics/momics/output/predictions/Manon_may26/Pombe_modulable_Manon_may26_no_maxpool",
    use_species_subdir = FALSE    # flat folder — same files for all species
  ),
  single_species = list(
    path             = "/pasteur/helix/scratch/adecugis/all_momics/momics/output/predictions/Manon_may26",
    use_species_subdir = TRUE     # has species subfolders
  )
)

GT_BASE_DIR  <- "/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species"
TSS_FOLDER   <- "/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/reference_genomes/5_yeasts/TSS_files_reference"
OUTDIR_BASE  <- "/pasteur/appa/scratch/adecugis/all_momics/momics/output/R/concatenated_profiles_Manon_May26"

KNOWN_TARGETS <- c("H3", "MNASE_15", "MNASE_25", "MNASE_35")

SPECIES_TSS <- list(
  C_glabrata   = file.path(TSS_FOLDER, "C_glabrata_TSS_forward.bed"),
  S_cerevisiae = file.path(TSS_FOLDER, "S_cerevisiae_TSS_forward.bed"),
  S_eubayanus  = file.path(TSS_FOLDER, "S_eubayanus_TSS_forward.bed"),
  S_paradoxus  = file.path(TSS_FOLDER, "S_paradoxus_TSS_forward.bed"),
  S_pombe      = file.path(TSS_FOLDER, "972hminus_TSS_forward.bed")
)

# =========================================================
# ANALYSIS FUNCTION
# =========================================================
run_analysis <- function(track_files, tss_file, species, target, outdir) {

  if (!file.exists(tss_file)) {
    message("❌ TSS file not found: ", tss_file); return(NULL) }

  features_gr <- tryCatch(import(tss_file),
    error = function(e) { message("❌ Failed to import TSS: ", e$message); NULL })
  if (is.null(features_gr) || length(features_gr) == 0) {
    message("❌ TSS file empty: ", tss_file); return(NULL) }

  bw_chroms  <- tryCatch(seqlevels(BigWigFile(unname(track_files[1]))),
                         error = function(e) character(0))
  overlap    <- intersect(bw_chroms, seqlevels(features_gr))
  if (length(overlap) == 0) {
    message("❌ No chrom overlap for ", target)
    message("   BW : ", paste(head(bw_chroms, 3),  collapse = ", "))
    message("   TSS: ", paste(head(seqlevels(features_gr), 3), collapse = ", "))
    return(NULL)
  }

  tracks <- BigWigFileList(setNames(as.list(track_files), names(track_files)))
  CE <- tryCatch(
    CoverageExperiment(tracks, list(TSSs = features_gr), width = kb_around,
                       ignore.strand = FALSE, BPPARAM = BiocParallel::SerialParam()),
    error = function(e) { message("❌ CoverageExperiment failed: ", e$message); NULL })
  if (is.null(CE)) return(NULL)

  agg <- aggregate(CE |> filter(features == "TSSs"))
  if (nrow(agg) == 0) { message("❌ Empty aggregation for ", target); return(NULL) }

  m <- as.data.frame(do.call(cbind, lapply(as.data.frame(assay(agg, "mean")), unlist)))

  # correlations — every prediction vs ground truth
  corr_tbl <- data.frame(Reference=character(), Track=character(),
                          Model=character(), Pearson=numeric(), Spearman=numeric(),
                          stringsAsFactors=FALSE)
  for (tname in names(track_files)) {
    if (grepl("^prediction_", tname)) {
      
      # ✅ Use known model names explicitly
      model_names_pattern <- paste(names(PRED_DIRS), collapse = "|")
      gt_name <- sub(
        paste0("^prediction_(", model_names_pattern, ")_"),
        "ground_truth_",
        tname
      )

      if (gt_name %in% names(track_files)) {
        x <- suppressWarnings(as.numeric(m[[gt_name]]))
        y <- suppressWarnings(as.numeric(m[[tname]]))
        model_label <- sub(paste0("^prediction_([^_]+)_.*"), "\\1", tname)
        corr_tbl <- rbind(corr_tbl, data.frame(
          Reference = gt_name, Track = tname, Model = model_label,
          Pearson  = round(suppressWarnings(cor(x, y, use="pairwise.complete.obs", method="pearson")),  3),
          Spearman = round(suppressWarnings(cor(x, y, use="pairwise.complete.obs", method="spearman")), 3),
          stringsAsFactors = FALSE))
      }
    }
  }
  print(corr_tbl)

  # plot — all tracks in one panel per target
  p <- ggplot(agg, aes(col = track)) +
    geom_aggrcoverage(linewidth = 0.6, alpha = 0.6) +
    facet_grid(track ~ ., scales = "fixed", switch = "y") +
    scale_color_manual(values = get_track_colors(length(unique(agg$track)))) +
    labs(title    = paste(species, "-", target),
         subtitle = paste("aggregated coverage (±", kb_around / 1000, "kb)"),
         x = "Distance from TSS (bp)", y = "Normalised signal", color = "Track") +
    theme_bw(base_size = 18) +
    theme(plot.title       = element_text(face = "bold", size = 18, hjust = 0.5),
          strip.text.y     = element_blank(),
          legend.position  = "top",
          panel.grid.minor = element_blank(),
          panel.spacing.y  = unit(0.4, "lines"))

  png_out <- file.path(outdir, paste0("coverage_", species, "_", target, "_", tag, "_TSSs.png"))
  ggsave(png_out, p, width = 8, height = 6, dpi = 300)
  message("  ✅ Plot saved: ", png_out)

  csv_out <- file.path(outdir, paste0("correlations_", species, "_", target, "_", tag, ".csv"))
  write.csv(corr_tbl, csv_out, row.names = FALSE)

  agg$Target      <- target
  corr_tbl$Target <- if (nrow(corr_tbl) > 0) target else character(0)
  return(list(agg = agg, corr_tbl = corr_tbl))
}

# =========================================================
# HELPER — extract target name from filename
# =========================================================
extract_target <- function(pred_file, known_targets) {
  stem   <- tools::file_path_sans_ext(pred_file)
  target <- sub("^(.+)_\\1$", "\\1", stem)
  if (target == stem) {
    matched <- known_targets[sapply(known_targets, function(t) grepl(t, stem, fixed=TRUE))]
    if (length(matched) == 0) return(NULL)
    target <- matched[which.max(nchar(matched))]
  }
  return(target)
}

# Add this helper near the top with the other helpers
find_species_subdir <- function(base_dir, species) {
  # try exact match
  exact <- file.path(base_dir, species)
  if (dir.exists(exact)) return(exact)

  # try case-insensitive match against all subdirs
  subdirs <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
  match   <- subdirs[tolower(subdirs) == tolower(species)]
  if (length(match) > 0) return(file.path(base_dir, match[1]))

  # try partial match (e.g. "glabrata" inside "c_glabrata")
  species_core <- tolower(gsub(".*_", "", species))
  match <- subdirs[grepl(species_core, tolower(subdirs), fixed = TRUE)]
  if (length(match) > 0) return(file.path(base_dir, match[1]))

  return(NULL)
}

# =========================================================
# MAIN LOOP — species → targets
# =========================================================
for (species in names(SPECIES_TSS)) {
  message("\n", paste(rep("=", 60), collapse=""))
  message("Species: ", species)

  gt_dir   <- file.path(GT_BASE_DIR, species)
  tss_file <- SPECIES_TSS[[species]]
  outdir   <- file.path(OUTDIR_BASE, species)

  if (!dir.exists(gt_dir)) { message("⚠️  No GT folder: ", gt_dir); next }
  dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

  gt_files <- list.files(gt_dir, pattern="\\.bw$", full.names=FALSE)
  if (length(gt_files) == 0) { message("⚠️  No GT files"); next }

  # collect pred files from BOTH model dirs for this species
  pred_files_by_model <- lapply(names(PRED_DIRS), function(model_name) {
      cfg <- PRED_DIRS[[model_name]]

      pred_dir <- if (cfg$use_species_subdir) {
        find_species_subdir(cfg$path, species)
      } else {
        cfg$path   # use flat folder directly
      }

      if (is.null(pred_dir) || !dir.exists(pred_dir)) {
        message("⚠️  No pred folder for model '", model_name,
                "' species '", species, "'")
        return(NULL)
      }

      message("  📁 Model '", model_name, "' pred dir: ", pred_dir)
      files <- list.files(pred_dir, pattern = "\\.bw$", full.names = FALSE)
      if (length(files) == 0) return(NULL)
      data.frame(model    = model_name,
                pred_dir = pred_dir,
                file     = files,
                stringsAsFactors = FALSE)
    })

  pred_df <- do.call(rbind, pred_files_by_model)

  if (is.null(pred_df) || nrow(pred_df) == 0) {
    message("⚠️  No prediction files found for any model — skipping ", species); next }

  # find all unique targets across both models
  pred_df$target <- sapply(pred_df$file, extract_target, known_targets=KNOWN_TARGETS)
  pred_df        <- pred_df[!sapply(pred_df$target, is.null), ]
  pred_df$target <- unlist(pred_df$target)
  all_targets    <- unique(pred_df$target)

  message("Targets found: ", paste(all_targets, collapse=", "))

  all_agg_list <- list()
  all_corr_tbl <- data.frame(Reference=character(), Track=character(), Model=character(),
                              Pearson=numeric(), Spearman=numeric(), Target=character(),
                              stringsAsFactors=FALSE)

  for (target in all_targets) {
    message("\n  Target: ", target)

    track_files <- c()

    # ground truth
    gt_matches <- grep(target, gt_files, value=TRUE, fixed=TRUE)
    if (length(gt_matches) > 0) {
      track_files[paste0("ground_truth_", target)] <- file.path(gt_dir, gt_matches[1])
      message("  ✅ GT: ", gt_matches[1])
    } else {
      message("  ⚠️  No GT for: ", target)
    }

    # one prediction entry per model
    for (model_name in names(PRED_DIRS)) {
      rows <- pred_df[pred_df$model == model_name & pred_df$target == target, ]
      if (nrow(rows) == 0) {
        message("  ⚠️  No prediction for model '", model_name, "' target '", target, "'")
        next
      }
      track_key <- paste0("prediction_", model_name, "_", target)
      track_files[track_key] <- file.path(rows$pred_dir[1], rows$file[1])
      message("  ✅ Pred (", model_name, "): ", rows$file[1])
    }

    missing <- !file.exists(unname(track_files))
    if (any(missing)) {
      message("  ❌ Missing: ", paste(track_files[missing], collapse=", ")); next }

    res <- run_analysis(track_files, tss_file, species, target, outdir)

    if (!is.null(res)) {
      all_agg_list[[target]] <- res$agg
      all_corr_tbl           <- rbind(all_corr_tbl, res$corr_tbl)
    }
  }

  # combined plot
  if (length(all_agg_list) > 0) {
    all_cols   <- unique(unlist(lapply(all_agg_list, function(a) names(as.data.frame(as_tibble(a))))))
    all_agg_df <- do.call(rbind, lapply(names(all_agg_list), function(t) {
      df <- as.data.frame(as_tibble(all_agg_list[[t]]))
      for (m in setdiff(all_cols, names(df))) df[[m]] <- NA
      df[all_cols]
    }))

    corr_labels <- all_corr_tbl %>%
      filter(!is.na(Pearson)) %>%
      group_by(Target) %>%
      summarise(
        label = paste(paste0(Model, ": r=", Pearson, " ρ=", Spearman), collapse="\n"),
        .groups = "drop"
      )

    p_all <- ggplot(all_agg_df, aes(col=track, group=track)) +
      geom_aggrcoverage(linewidth=0.6, alpha=0.6) +
      geom_text(data=corr_labels, aes(label=label),
                x=-Inf, y=Inf, hjust=-0.05, vjust=1.3,
                size=7, color="grey20", inherit.aes=FALSE) +
      facet_wrap(~Target, ncol=4, scales="fixed") +
      scale_color_manual(values = get_track_colors(length(unique(all_agg_df$track)))) +
      labs(title    = paste(species, "— all targets —", tag),
           subtitle = paste("aggregated coverage (±", kb_around/1000, "kb)"),
           x="Distance from TSS (bp)", y="Normalised signal", color="Track") +
      theme_bw(base_size=18) +
      theme(strip.text=element_text(size=18, face="bold"),
            legend.position="top", panel.grid.minor=element_blank(), legend.text=element_text(size=17))  

    ggsave(file.path(outdir, paste0("coverage_", species, "_", tag, "_ALL.png")),
           p_all, width=20, height=12, dpi=300)
    write.csv(all_corr_tbl,
              file.path(outdir, paste0("correlations_", species, "_", tag, "_ALL.csv")),
              row.names=FALSE)
    message("✅ Combined outputs saved for ", species)
  }
}

message("\n✅ All species complete.")