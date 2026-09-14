# ==============================================================================
# RGE AGGREGATION — combine leave-one-variable-out predictions with the
# full-model predictions to compute the RGE (Rank Graduation Explainability)
# score for every (dataset, variable, model) triple.
#
#   RGE = 1 - RGEstar(yhat_full, yhat_minus_variable)
#
# RGEstar is the same Rank Graduation Ratio (RGR)-style statistic used
# elsewhere in the paper (fig:bench_rgr): it compares two prediction vectors
# via their induced rankings. RGE close to 1 means removing the variable
# changes the model's ranking of predictions a lot (the model relies heavily
# on that variable); RGE close to 0 means removing it barely matters.
#
# Inputs:
#   pred_RGE_benchmark/                    (this script's sibling folder --
#                                            produced by run_rge_predictions.R)
#   ../vettori_predizioni_benchmark_ADAC/  (full-model rep-1 predictions --
#                                            produced by the run_benchmark_*.R
#                                            scripts in Benchmark/)
# Outputs:
#   output_RGE/df_RGE_long.rds
#   output_RGE/df_RGE_wide.rds   (this is what fig6_rge_explainability.R reads)
#
# NOTE: "from scratch" aggregation script (slow: one file read per
# dataset x variable x model). NOT required to reproduce fig:rge_vs_sword --
# fig6_rge_explainability.R already reads the deposited output_RGE/df_RGE_wide.rds
# directly. Checkpointed by chunk so a long run can be resumed.
# ==============================================================================

USE_PARALLEL <- TRUE
CHUNK_SIZE   <- 500

library(dplyr)
library(tidyr)

## ---- Path-robust header (Rscript and RStudio) ----
.get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg)) return(dirname(normalizePath(file_arg)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}
SCRIPT_DIR <- .get_script_dir()                              # .../analisi_RGE_explainability
BENCH_DIR  <- normalizePath(file.path(SCRIPT_DIR, ".."))     # .../Benchmark

pred_rge_dir    <- file.path(SCRIPT_DIR, "pred_RGE_benchmark")
pred_full_dir   <- file.path(BENCH_DIR, "vettori_predizioni_benchmark_ADAC")
out_dir         <- file.path(SCRIPT_DIR, "output_RGE")
checkpoint_file <- file.path(out_dir, "checkpoint_rge_vals.rds")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

N_CORES <- max(1L, parallel::detectCores() - 5)

MODELS <- c("SWORD", "RF", "ODRF", "ROT", "EXTree", "RRF",
            "aorsf", "aorsfNET", "SPORF", "CF", "SVMlin")

# ==============================================================================
# 1. RGEstar (same Rank Graduation Ratio statistic as fig:bench_rgr)
# ==============================================================================

RGEstar <- function(yhat, yhat_pert) {
  set.seed(1)
  ryhat_pert <- rank(round(yhat_pert, 4), ties.method = "min")
  support <- tapply(yhat, ryhat_pert, mean)
  rord <- c(1:length(yhat))
  for (jj in 1:length(yhat)) {
    rord[jj] <- support[names(support) == ryhat_pert[jj]]
  }
  yhatstar <- rord[order(yhat_pert)]
  I    <- 1:length(yhatstar)
  conc <- 2 * sum(I * yhatstar)
  dec  <- 2 * sum(I * sort(yhat, decreasing = TRUE))
  inc  <- 2 * sum(I * sort(yhat))
  (conc - dec) / (inc - dec)
}

# ==============================================================================
# 2. Discover all pred_RGE files and build the task index
# ==============================================================================

all_files <- list.files(pred_rge_dir, pattern = "_VARIAB_.*_rep1\\.rds$", full.names = FALSE)
cat("pred_RGE files found:", length(all_files), "\n")

parse_rge_fname <- function(fname) {
  core <- sub("^pred_RGE_", "", sub("_rep1\\.rds$", "", fname))
  for (m in MODELS) {
    prefix <- paste0(m, "_")
    if (startsWith(core, prefix)) {
      rest  <- sub(paste0("^", m, "_"), "", core)
      parts <- strsplit(rest, "_VARIAB_", fixed = TRUE)[[1]]
      if (length(parts) != 2) return(NULL)
      return(data.frame(model = m, safe_ds = parts[1], var_name = parts[2],
                        fname = fname, stringsAsFactors = FALSE))
    }
  }
  NULL
}

index_df <- do.call(rbind, lapply(all_files, parse_rge_fname))
index_df <- index_df[!is.na(index_df$model), ]
n_rows   <- nrow(index_df)
cat("Valid index rows:", n_rows, "\n")
cat("Unique datasets:", dplyr::n_distinct(index_df$safe_ds), "\n")
cat("Models:", paste(sort(unique(index_df$model)), collapse = ", "), "\n\n")

# ==============================================================================
# 3. Pre-load all full-model predictions into memory
# ==============================================================================

cat("Pre-loading full-model predictions...\n")
unique_keys <- unique(paste0(index_df$model, "__", index_df$safe_ds))
full_preds  <- vector("list", length(unique_keys))
names(full_preds) <- unique_keys

for (key in unique_keys) {
  parts      <- strsplit(key, "__", fixed = TRUE)[[1]]
  m          <- parts[1]
  safe       <- parts[2]
  model_full <- if (m == "SWORD") "SWORDSCALE" else m   # SWORD's full-model files use the "SWORDSCALE" label
  f <- file.path(pred_full_dir, paste0("pred_", model_full, "_", safe, "_rep1.rds"))
  full_preds[[key]] <- if (file.exists(f)) {
    tryCatch(as.numeric(readRDS(f)), error = function(e) NULL)
  } else NULL
}

n_loaded <- sum(!sapply(full_preds, is.null))
cat("Full-model predictions loaded:", n_loaded, "/", length(unique_keys), "\n\n")

# ==============================================================================
# 4. Load checkpoint if present (resume where it left off)
# ==============================================================================

results_df <- data.frame(
  model = index_df$model, dataset = index_df$safe_ds,
  variabile = index_df$var_name, RGE = NA_real_,
  stringsAsFactors = FALSE
)

if (file.exists(checkpoint_file)) {
  results_df <- readRDS(checkpoint_file)
  cat(sprintf("Checkpoint found: %d / %d values already computed. Resuming...\n\n",
              sum(!is.na(results_df$RGE)), n_rows))
} else {
  cat("No checkpoint found. Starting from scratch.\n\n")
}

todo <- which(is.na(results_df$RGE))
cat(sprintf("Indices left to compute: %d\n\n", length(todo)))

# ==============================================================================
# 5. Compute RGE, chunk by chunk (checkpointed)
# ==============================================================================

chunks   <- split(todo, ceiling(seq_along(todo) / CHUNK_SIZE))
n_chunks <- length(chunks)

if (USE_PARALLEL) {
  library(parallel); library(foreach); library(doParallel)
  cat(sprintf("Parallel: %d cores\n\n", N_CORES))
}

compute_one <- function(i) {
  m       <- as.character(index_df$model[i])
  safe_ds <- as.character(index_df$safe_ds[i])
  var_nm  <- as.character(index_df$var_name[i])

  yhat_minus <- tryCatch(as.numeric(readRDS(file.path(pred_rge_dir, index_df$fname[i]))),
                        error = function(e) NULL)
  rge_val <- NA_real_
  if (!is.null(yhat_minus) && length(yhat_minus) > 0) {
    yhat_full <- full_preds[[paste0(m, "__", safe_ds)]]
    if (!is.null(yhat_full) && length(yhat_full) > 0) {
      rge_val <- tryCatch(1 - RGEstar(yhat_full, yhat_minus), error = function(e) NA_real_)
    }
  }
  data.frame(model = m, dataset = safe_ds, variabile = var_nm, RGE = rge_val, stringsAsFactors = FALSE)
}

for (ch in seq_len(n_chunks)) {
  idx_chunk <- chunks[[ch]]
  cat(sprintf("--- Chunk %d / %d (indices %d-%d, %d rows) ---\n",
              ch, n_chunks, idx_chunk[1], idx_chunk[length(idx_chunk)], length(idx_chunk)))

  if (USE_PARALLEL) {
    cl <- makeCluster(N_CORES)
    registerDoParallel(cl)
    clusterExport(cl, varlist = c("pred_rge_dir", "full_preds", "RGEstar", "index_df", "compute_one"),
                  envir = environment())
    chunk_df <- foreach(i = idx_chunk, .combine = rbind, .packages = character(0)) %dopar% compute_one(i)
    stopCluster(cl)
  } else {
    chunk_df <- do.call(rbind, lapply(idx_chunk, compute_one))
  }

  results_df$RGE[idx_chunk] <- chunk_df$RGE
  saveRDS(results_df, checkpoint_file)
  cat(sprintf("    Checkpoint saved: %d / %d total (NA: %d)\n\n",
              sum(!is.na(results_df$RGE)), n_rows, sum(is.na(results_df$RGE))))
}

cat(sprintf("Done. RGE computed: %d | NA: %d\n\n",
            sum(!is.na(results_df$RGE)), sum(is.na(results_df$RGE))))

# ==============================================================================
# 6. Build long and wide data frames from the computed RGE values
# ==============================================================================

df_long <- results_df %>%
  dplyr::select(dataset = dataset, var_name = variabile, model, RGE) %>%
  arrange(dataset, var_name, model)

df_wide <- df_long %>%
  pivot_wider(names_from = model, values_from = RGE, names_glue = "RGE_{model}") %>%
  arrange(dataset, var_name)

cat("df_wide dimensions:", nrow(df_wide), "rows x", ncol(df_wide), "columns\n\n")

# ==============================================================================
# 7. Save
# ==============================================================================

saveRDS(df_long, file.path(out_dir, "df_RGE_long.rds"))
saveRDS(df_wide, file.path(out_dir, "df_RGE_wide.rds"))
cat("Saved: output_RGE/df_RGE_long.rds and df_RGE_wide.rds\n")
cat("\nEND aggregate_rge_scores.R\n")
