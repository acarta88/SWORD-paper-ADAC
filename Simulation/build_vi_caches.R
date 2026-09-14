# ==============================================================================
# Build light VI caches (seed_1 only) from the raw per-dataset simulation files
#
# Purpose:
#   `run_simulation_scale.R` and `run_simulation_RF.R` already compute variable
#   importance (OIW-VI for SWORD, MDI for RF) as part of their normal run, but
#   only for seed_1, and they store it inside each per-dataset raw .rds file
#   (field $VI). This script scans those raw files, keeps only the recommended
#   SWORD config (W_scheme=scale | Pearson | rand_ntopcor=TRUE | Cost_C=1 |
#   rf_var_frac=100%), and compacts everything into two small cache files:
#     - risultati_full_simulation_scale/VI_SWORD_seed1_scale_cached.rds
#     - risultati_RF_simulation/VI_RF_seed1_cached.rds
#   These are exactly the caches consumed by fig1_variable_importance.R (which contains the
#   same extraction logic inline as a fallback, in case the caches are missing
#   but the raw files are still present).
#
# When to run this script:
#   Only needed if you regenerate the raw per-dataset files from scratch via
#   run_simulation_scale.R / run_simulation_RF.R and want to refresh the two
#   light caches shipped in this repo. The raw per-dataset .rds files are NOT
#   included in this repo (only the two caches are, to keep it lightweight) —
#   running this script as-is on a fresh clone will find 0 files and exit
#   without overwriting the shipped caches.
# ==============================================================================

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
SCALE_DIR  <- file.path(SCRIPT_DIR, "risultati_full_simulation_scale")
RF_DIR     <- file.path(SCRIPT_DIR, "risultati_RF_simulation")


# ==============================================================================
# 1. VI SWORD (scale, seed_1, recommended combo only)
#    File name pattern: scale__seed_1__<n_obs>_n_obs__<feat>_n_feat__...
#    val: [NA(scale), seed, n_obs, features, noise_prop, nonlin, cat_prop, error_sd]
# ==============================================================================

file_rds_sword <- list.files(SCALE_DIR, pattern = "\\.rds$", full.names = TRUE)
file_rds_sword <- file_rds_sword[basename(file_rds_sword) != "risultati_completi.rds"]
file_rds_sword <- file_rds_sword[!grepl("^VI_.*_cached\\.rds$", basename(file_rds_sword))]
file_rds_sword <- file_rds_sword[grepl("^scale__seed_1__", basename(file_rds_sword))]
cat("SWORD scale seed_1 raw files found:", length(file_rds_sword), "\n")

if (length(file_rds_sword) == 0) {

  cat("No raw SWORD files found - nothing to rebuild.\n")
  cat("(Re-run run_simulation_scale.R first if you want to regenerate the cache.)\n")

} else {

  vi_sword_raw <- do.call(rbind, Filter(Negate(is.null), lapply(file_rds_sword, function(f) {
    nm     <- gsub("\\.rds$", "", basename(f))
    el     <- strsplit(nm, "__")[[1]]
    val    <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", el)))
    n_feat <- val[4]

    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(obj)) return(NULL)

    combo_list <- if (!is.null(obj$combos)) obj$combos else
      Filter(function(x) is.list(x) && !is.null(x$res), obj)
    if (length(combo_list) == 0) return(NULL)

    do.call(rbind, Filter(Negate(is.null), lapply(combo_list, function(x) {
      rel   <- x$res$relation
      rnd   <- as.logical(x$res$rand_ntopcor)
      costC <- as.numeric(as.character(x$res$Cost_C))
      frac  <- x$res$rf_var / n_feat

      # Keep only the recommended config: Pearson | Random gamma | C=1 | alpha=100%
      if (!isTRUE(rel == "Pearson")) return(NULL)
      if (!isTRUE(rnd == TRUE))     return(NULL)
      if (!isTRUE(costC == 1))      return(NULL)
      if (frac < 0.75)               return(NULL)

      vi <- x$VI
      if (is.null(vi) || length(vi) == 0 || is.null(names(vi))) return(NULL)

      data.frame(
        variable   = names(vi),
        VI_SWORD   = as.numeric(vi),
        n_obs      = val[3],
        features   = n_feat,
        noise_prop = val[5],
        nonlin     = val[6],
        cat_prop   = val[7],
        error_sd   = val[8],
        stringsAsFactors = FALSE
      )
    })))
  })))

  cat("  VI SWORD rows extracted:", nrow(vi_sword_raw), "\n")
  saveRDS(vi_sword_raw, file.path(SCALE_DIR, "VI_SWORD_seed1_scale_cached.rds"))
  cat("  Saved:", file.path(SCALE_DIR, "VI_SWORD_seed1_scale_cached.rds"), "\n")
}


# ==============================================================================
# 2. VI RF (seed_1 only)
#    File name pattern: seed_1__<n_obs>_n_obs__<feat>_n_feat__...
#    val: [seed, n_obs, features, noise_prop, nonlin, cat_prop, error_sd]
# ==============================================================================

file_rds_rf <- list.files(RF_DIR, pattern = "\\.rds$", full.names = TRUE)
file_rds_rf <- file_rds_rf[basename(file_rds_rf) != "risultati_completi.rds"]
file_rds_rf <- file_rds_rf[!grepl("^VI_.*_cached\\.rds$", basename(file_rds_rf))]
file_rds_rf <- file_rds_rf[grepl("^seed_1__", basename(file_rds_rf))]
cat("\nRF seed_1 raw files found:", length(file_rds_rf), "\n")

if (length(file_rds_rf) == 0) {

  cat("No raw RF files found - nothing to rebuild.\n")
  cat("(Re-run run_simulation_RF.R first if you want to regenerate the cache.)\n")

} else {

  vi_rf_raw <- do.call(rbind, Filter(Negate(is.null), lapply(file_rds_rf, function(f) {
    nm  <- gsub("\\.rds$", "", basename(f))
    el  <- strsplit(nm, "__")[[1]]
    val <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", el)))

    res <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(res)) return(NULL)

    do.call(rbind, Filter(Negate(is.null), lapply(res, function(x) {
      vi <- x$VI
      if (is.null(vi)) return(NULL)
      if (is.matrix(vi) || is.data.frame(vi)) {
        vi_names <- rownames(vi)
        vi_vals  <- as.numeric(vi[, 1])
      } else {
        vi_names <- names(vi)
        vi_vals  <- as.numeric(vi)
      }
      if (is.null(vi_names) || length(vi_names) == 0) return(NULL)
      # RF MDI may include small negative values -> clip to 0
      vi_vals <- pmax(vi_vals, 0)

      data.frame(
        variable   = vi_names,
        VI_RF      = vi_vals,
        n_obs      = val[2],
        features   = val[3],
        noise_prop = val[4],
        nonlin     = val[5],
        cat_prop   = val[6],
        error_sd   = val[7],
        stringsAsFactors = FALSE
      )
    })))
  })))

  cat("  VI RF rows extracted:", nrow(vi_rf_raw), "\n")
  saveRDS(vi_rf_raw, file.path(RF_DIR, "VI_RF_seed1_cached.rds"))
  cat("  Saved:", file.path(RF_DIR, "VI_RF_seed1_cached.rds"), "\n")
}

cat("\nDone.\n")
