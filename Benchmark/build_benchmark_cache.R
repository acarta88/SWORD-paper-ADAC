# ==============================================================================
# Build the light benchmark cache from the raw per-dataset/per-model files
#
# Purpose:
#   fig3_fig4_benchmark_rmse_rga.R needs one merged table (df_all: SWORDSCALE + the
#   10 competitor models, harmonised to a common set of columns) to produce
#   fig:bench_res_all and fig:bench_res_RGA. This script performs exactly that
#   load-and-merge step (same logic as the "cache not found" branch inside
#   fig3_fig4_benchmark_rmse_rga.R) and saves the result as a small standalone cache:
#     - results_cache/df_benchmark248_cached.rds
#
# When to run this script:
#   Only needed if you regenerate the raw per-dataset files from scratch via
#   run_benchmark_SWORD.R / run_benchmark_{RF,ODRF,ROT,EXTree,RRF,aorsf,
#   aorsfNET,SPORF,CF,SVMlin}.R and want to refresh the cache shipped in this
#   repo. The raw files themselves (11362 of them, ~31 MB) are NOT included in
#   this repo, only the cache is — running this script as-is on a fresh clone
#   will find 0 files and exit without overwriting the shipped cache.
# ==============================================================================

library(dplyr)

## ---- Path-robust header (Rscript and RStudio) ----
.get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg)) return(dirname(normalizePath(file_arg)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}
SCRIPT_DIR <- .get_script_dir()
PROJ    <- normalizePath(file.path(SCRIPT_DIR, ".."))
CACHE   <- file.path(PROJ, "results_cache")
dir_res <- file.path(SCRIPT_DIR, "risultati_finali_benchmark_ADAC")

benchmark_cache_path <- file.path(CACHE, "df_benchmark248_cached.rds")


# ==============================================================================
# 1. LOAD SWORDSCALE RESULTS
#    File batch:  risultati_ADAC_SWORDSCALE_<ds>_1_50
# ==============================================================================

swordscale_files_batch <- list.files(dir_res,
                                pattern = "^risultati_ADAC_SWORDSCALE_.*_[0-9]+_[0-9]+\\.rds$",
                                full.names = TRUE)
swordscale_files_rep <- list.files(dir_res,
                              pattern = "^risultati_ADAC_SWORDSCALE_.*_rep[0-9]+\\.rds$",
                              full.names = TRUE)
swordscale_files_batch <- swordscale_files_batch[!grepl("_ALL_", swordscale_files_batch)]

n_raw <- length(swordscale_files_batch) + length(swordscale_files_rep)
cat("SWORDSCALE raw files found:", n_raw, "\n")

if (n_raw == 0) {

  cat("No raw benchmark files found - nothing to rebuild.\n")
  cat("(Re-run run_benchmark_SWORD.R and the 10 competitor scripts first if you want to regenerate the cache.)\n")

} else {

  cat("  File batch:", length(swordscale_files_batch), "\n")
  cat("  Single-rep files:", length(swordscale_files_rep), "\n")

  read_safe <- function(f) {
    tryCatch(readRDS(f), error = function(e) {
      warning("Could not read: ", f); NULL
    })
  }

  list_swordscale_batch <- lapply(swordscale_files_batch, read_safe)
  list_swordscale_rep   <- lapply(swordscale_files_rep,   read_safe)

  df_swordscale_raw <- do.call(rbind, c(list_swordscale_batch, list_swordscale_rep))
  df_swordscale_raw$model <- "SWORD"
  df_swordscale_raw$spec  <- "SWORD"

  # Deduplication: keep a single row per (dataset, repetition)
  df_swordscale <- df_swordscale_raw %>%
    arrange(dataset, repetition) %>%
    group_by(dataset, repetition) %>%
    slice_tail(n = 1) %>%
    ungroup()

  cat("  SWORDSCALE rows after dedup:", nrow(df_swordscale), "\n")
  cat("  Dataset SWORDSCALE:", n_distinct(df_swordscale$dataset), "\n\n")


  # ============================================================================
  # 2. LOAD OTHER MODELS' RESULTS
  # ============================================================================
  cat("=== Loading other models ===\n")

  OTHER_MODELS <- c("RF", "ODRF", "ROT", "EXTree", "RRF",
                    "aorsfNET", "aorsf", "SPORF", "CF", "SVMlin")

  load_model <- function(model_name) {
    pattern <- paste0("^risultati_", model_name, "_ADAC_.*\\.rds$")
    files   <- list.files(dir_res, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) { cat("  No file for:", model_name, "\n"); return(NULL) }
    files <- files[!grepl("_ALL_", files)]
    lst <- lapply(files, read_safe)
    df  <- do.call(rbind, lst)
    if (!is.null(df) && nrow(df) > 0) {
      df <- df[!duplicated(df[, c("dataset", "repetition")]), ]
    }
    if (!"model" %in% names(df)) df$model <- model_name
    if (!"spec"  %in% names(df)) df$spec  <- model_name
    df$model <- model_name
    df$spec  <- model_name
    cat("  ", model_name, ":", nrow(df), "rows,",
        n_distinct(df$dataset), "dataset\n")
    df
  }

  list_other <- lapply(OTHER_MODELS, load_model)
  names(list_other) <- OTHER_MODELS


  # ============================================================================
  # 3. MERGE INTO A SINGLE DATA FRAME AND SAVE
  # ============================================================================
  cat("\n=== Merge ===\n")

  CORE_COLS <- c("dataset", "repetition", "RMSE", "MAE", "R2", "RGA",
                 "sd_y_train", "sd_y_test", "ncolumn",
                 "cor_tree", "strength_tree", "time_forest",
                 "model", "spec", "ntree")

  harmonise <- function(df, cols = CORE_COLS) {
    if (is.null(df)) return(NULL)
    for (col in cols) {
      if (!col %in% names(df)) df[[col]] <- NA
    }
    df[, cols, drop = FALSE]
  }

  df_all <- do.call(rbind,
                    c(list(harmonise(df_swordscale)),
                      lapply(list_other, harmonise)))

  cat("  Total rows:", nrow(df_all), "\n")
  cat("  Total datasets:", n_distinct(df_all$dataset), "\n")
  cat("  Models:", paste(unique(df_all$model), collapse = ", "), "\n\n")

  if (!dir.exists(CACHE)) dir.create(CACHE, recursive = TRUE)
  saveRDS(df_all, benchmark_cache_path)
  cat("Saved:", benchmark_cache_path, "\n")
}
