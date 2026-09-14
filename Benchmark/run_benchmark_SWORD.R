# ==============================================================================
# FULL BENCHMARK GENERATION — SWORD (scale scheme) on the 248 real datasets
# Standard cycle (50 repetitions) + perturbed-variables robustness cycle
#
# NOTE: this is the "from scratch" generation script (slow: hours, meant for a
# many-core machine). It is included for full transparency/provenance and is
# NOT required to reproduce any paper figure/table (fig3_fig4_benchmark_rmse_rga.R,
# fig5_benchmark_rgr.R, and the RGE scripts read the cached/deposited results
# in results_cache/). Run this only to regenerate the SWORD benchmark from the
# raw 248 real datasets.
#
# Recipe:
#   nmin=5, minleaf=2, cp=0, n_perc=1, n_topCor=2, threshold_COR=1, cost_C=1,
#   m=100, rf_var=ncol(x_train) (i.e. 100% of the retained features),
#   Weight_Scheme="scale", rand_ntopcor=TRUE, relation="Pearson" (default),
#   OOB=FALSE. Train/test split: 70/30 via caret::createDataPartition with
#   set.seed(10), 200 pre-generated partitions (or 3500-row random subsamples
#   for datasets with > 5000 rows), repetitions 1:50 index into that matrix.
# ==============================================================================

library(fastDummies)
library(caret)
library(Metrics)
library(doParallel)
library(foreach)

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
PROJ       <- normalizePath(file.path(SCRIPT_DIR, ".."))
DATA       <- file.path(PROJ, "data")

source(file.path(PROJ, "R", "SVM.ROT_Functions_ADAC.R"))

load(file.path(DATA, "benchmark_248_dataset_list_ADAC.rda"))
names_248datasets <- names(benchmark_248_dataset_list)

# Output folders (created if missing; not shipped in this repo, see header above)
DIR_RESULTS  <- file.path(SCRIPT_DIR, "risultati_finali_benchmark_ADAC")
DIR_PRED     <- file.path(SCRIPT_DIR, "vettori_predizioni_benchmark_ADAC")
DIR_PERTURB  <- file.path(SCRIPT_DIR, "vettori_perturbed_benchmark_ADAC")
dir.create(DIR_RESULTS, showWarnings = FALSE)
dir.create(DIR_PRED,    showWarnings = FALSE)
dir.create(DIR_PERTURB, showWarnings = FALSE)

start_rep    <- 1
end_rep      <- 50
end_rep_pert <- 50

########## Helper functions ##########

sanitize_colnames <- function(nomi) {
  nomi <- gsub("[><=]", "_", nomi)
  nomi <- gsub("[^[:alnum:]_]", "_", nomi)
  nomi <- gsub("^[0-9]", "X", nomi)
  nomi <- gsub("weight", "WEIGHT", nomi)
  nomi <- gsub("time",   "TIME",   nomi)
  nomi <- make.names(nomi, unique = TRUE)
  return(nomi)
}

rename_except_target <- function(df, target_name = "target") {
  if (!(target_name %in% colnames(df))) {
    warning("Target column not found. Proceeding without exclusion.")
    target_name <- NULL
  }
  cols_to_rename <- setdiff(colnames(df), target_name)
  generate_names <- function(n) {
    base <- LETTERS
    if (n <= length(base)) return(paste0("var", base[1:n]))
    comb <- c(base, outer(base, base, paste0))
    return(paste0("var", comb[1:n]))
  }
  new_names <- generate_names(length(cols_to_rename))
  colnames(df)[colnames(df) %in% cols_to_rename] <- new_names
  return(df)
}

# Rank Graduation Accuracy (RGA), rank-based measure in [0, 1]:
#   1   = perfect concordance (predicted ranking = true ranking)
#   0.5 = ~random behaviour
#   0   = perfect discordance (predicted ranking = reverse of the true one)
RGA_metric <- function(y_true, y_pred) {
  y_true <- as.numeric(y_true)
  y_pred <- as.numeric(y_pred)
  ok     <- is.finite(y_true) & is.finite(y_pred)
  y_true <- y_true[ok]
  y_pred <- y_pred[ok]
  n <- length(y_true)
  if (n < 2) return(NA_real_)
  r_true     <- order(y_true)
  r_true_rev <- order(y_true, decreasing = TRUE)
  r_pred     <- order(y_pred)
  S_pred  <- sum((1:n) * y_true[r_pred])
  S_best  <- sum((1:n) * y_true[r_true])
  S_worst <- sum((1:n) * y_true[r_true_rev])
  num <- S_pred - S_worst
  den <- S_best - S_worst
  if (den == 0) return(NA_real_)
  num / den
}

# Deterministic perturbation (5%-95% swap), as in Babaei et al. 2025
perturb_variables <- function(x_to_perturb, perturbation_percentage = 0.05) {
  x_pert  <- as.data.frame(x_to_perturb)
  n_total <- nrow(x_pert)
  for (j in seq_len(ncol(x_pert))) {
    col <- as.numeric(x_pert[, j])
    if (length(unique(col)) <= 10) next   # skip categorical/low-cardinality
    sorted_idx <- order(col)
    p5_n    <- as.integer(ceiling(perturbation_percentage * n_total))
    p95_n   <- as.integer(ceiling((1 - perturbation_percentage) * n_total))
    n_upper <- n_total - p95_n
    n       <- min(p5_n, n_upper)
    if (n == 0) next
    lower_idx <- sorted_idx[seq_len(n)]
    upper_idx <- rev(sorted_idx[(n_total - n + 1):n_total])
    col_new            <- col
    col_new[lower_idx] <- col[upper_idx]
    col_new[upper_idx] <- col[lower_idx]
    x_pert[, j]        <- col_new
  }
  x_pert
}

########## Parallel backend ##########

n_workers <- min(10, parallel::detectCores() - 1L)
cat("Workers:", n_workers, "\n")

if (.Platform$OS.type == "windows") {

  cl <- makeCluster(n_workers)
  registerDoParallel(cl)

  clusterExport(cl, varlist = "PROJ")
  clusterEvalQ(cl, {
    Sys.setenv(OMP_NUM_THREADS      = "1")
    Sys.setenv(OPENBLAS_NUM_THREADS = "1")
    Sys.setenv(MKL_NUM_THREADS      = "1")
    source(file.path(PROJ, "R", "SVM.ROT_Functions_ADAC.R"))
    library(Metrics)
    library(fastDummies)
    library(caret)
  })
  clusterExport(cl, varlist = c(
    "benchmark_248_dataset_list", "DIR_RESULTS", "DIR_PRED", "DIR_PERTURB",
    "RGA_metric", "sanitize_colnames", "rename_except_target", "perturb_variables"
  ))

} else {

  library(doMC)
  Sys.setenv(OMP_NUM_THREADS      = "1")
  Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  Sys.setenv(MKL_NUM_THREADS      = "1")
  registerDoMC(cores = n_workers)

}

############################################################
# CYCLE 1 — SWORD standard
# foreach over all (dataset x rep) combinations
############################################################

combos_std <- expand.grid(
  ds  = names_248datasets,
  rep = seq(start_rep, end_rep),
  stringsAsFactors = FALSE
)
cat("\nCYCLE 1: SWORD standard --", nrow(combos_std), "combinations (ds x rep)\n")

risultati_raw <- foreach(
  ds_name = combos_std$ds,
  rep     = combos_std$rep,
  .options.multicore = list(preschedule = FALSE),
  .errorhandling = "remove"
) %dopar% {
  tryCatch({

    safe_ds  <- gsub("[^A-Za-z0-9_]", "_", ds_name)
    rep_file <- file.path(DIR_RESULTS,
                          paste0("risultati_ADAC_SWORDSCALE_", safe_ds,
                                 "_rep", rep, ".rds"))

    # -- Per-rep checkpoint --
    if (file.exists(rep_file)) {
      cat("Cached result:", ds_name, "rep", rep, "\n")
      return(readRDS(rep_file))
    }

    pred_file <- file.path(DIR_PRED, paste0("pred_SWORDSCALE_", safe_ds, "_rep", rep, ".rds"))
    true_file <- file.path(DIR_PRED, paste0("true_", safe_ds, "_rep", rep, ".rds"))

    # -- Preprocessing: deterministic, fast, always executed --
    df <- benchmark_248_dataset_list[[ds_name]]
    colnames(df) <- gsub("-", "_", colnames(df))
    colnames(df) <- gsub(" ", "_", colnames(df))
    if (any(sapply(df, is.factor) | sapply(df, is.character))) {
      df <- fastDummies::dummy_cols(df, remove_selected_columns = TRUE,
                                    remove_first_dummy = TRUE)
    }
    colnames(df) <- sanitize_colnames(colnames(df))
    df <- rename_except_target(df)

    y <- df$target
    x <- df[, setdiff(names(df), "target")]

    if (nrow(df) > 5000) {
      partition_matrix <- matrix(NA, nrow = 3500, ncol = 200)
      set.seed(10)
      for (ii in 1:200) partition_matrix[, ii] <- sample(nrow(df), 3500)
    } else {
      set.seed(10)
      partition_matrix <- caret::createDataPartition(df$target, p = 0.7,
                                                     times = 200, list = FALSE)
    }

    idx_train <- partition_matrix[, rep]
    x_train   <- x[idx_train,  , drop = FALSE]
    y_train   <- y[idx_train]
    x_test    <- x[-idx_train, , drop = FALSE]
    y_test    <- y[-idx_train]

    constant_cols <- sapply(x_train, function(col) length(unique(col)) == 1)
    x_train <- x_train[, !constant_cols, drop = FALSE]
    x_test  <- x_test[,  !constant_cols, drop = FALSE]
    rf_val  <- ncol(x_train)

    cat("Compute:", ds_name, "rep", rep, "\n")
    start_time <- Sys.time()
    model <- SVM.ROT.RF.OOB(
      Covariates    = x_train,
      y             = y_train,
      nmin          = 5,
      minleaf       = 2,
      cp            = 0,
      n_perc        = 1,
      n_topCor      = 2,
      threshold_COR = 1,
      cost_C        = 1,
      m             = 100,
      rf_var        = rf_val,
      Weight_Scheme = "scale",
      rand_ntopcor  = TRUE,
      parallel      = FALSE,
      OOB           = FALSE,
      verbose       = TRUE
    )
    end_time    <- Sys.time()
    time_forest <- as.numeric(difftime(end_time, start_time, units = "secs"))

    pred_matrix_all <- SVM.ROT.PRED.RF_all(x_test, model)
    pred_test       <- rowMeans(pred_matrix_all, na.rm = TRUE)

    cor_mat       <- cor(pred_matrix_all, use = "pairwise.complete.obs")
    cor_tree      <- mean(cor_mat[upper.tri(cor_mat)], na.rm = TRUE)
    r2_trees      <- apply(pred_matrix_all, 2, function(p)
      (cor(y_test, p, use = "pairwise.complete.obs"))^2)
    strength_tree  <- mean(r2_trees, na.rm = TRUE)
    time_tree_mean <- mean(model$time_tree, na.rm = TRUE)

    saveRDS(as.numeric(pred_test), pred_file)
    saveRDS(as.numeric(y_test),    true_file)

    # -- Metrics --
    rmse_val <- Metrics::rmse(y_test, pred_test)
    mae_val  <- Metrics::mae(y_test, pred_test)
    r2_val   <- (cor(y_test, pred_test))^2
    RGA_val  <- RGA_metric(y_test, pred_test)

    df_row <- data.frame(
      dataset        = ds_name,
      repetition     = rep,
      sd_y_train     = sd(y_train, na.rm = TRUE),
      sd_y_test      = sd(y_test,  na.rm = TRUE),
      ncolumn        = ncol(x),
      RMSE           = rmse_val,
      MAE            = mae_val,
      R2             = r2_val,
      RGA            = RGA_val,
      cor_tree       = cor_tree,
      strength_tree  = strength_tree,
      time_forest    = time_forest,
      time_tree_mean = time_tree_mean,
      mtry           = "SWORDSCALE",
      ntree          = 100
    )

    # -- Save per-rep checkpoint --
    saveRDS(df_row, rep_file)
    df_row

  }, error = function(e) {
    cat("Error (std):", ds_name, "rep", rep, ":", e$message, "\n")
    NULL
  })
}

# -- Aggregate the per-rep results into one bundled file per dataset --
cat("\nAggregating per-repetition results into one file per dataset...\n")
df_std_all <- do.call(rbind, risultati_raw)
for (ds_name in names_248datasets) {
  safe_ds <- gsub("[^A-Za-z0-9_]", "_", ds_name)
  df_ds   <- df_std_all[df_std_all$dataset == ds_name, ]
  if (nrow(df_ds) == 0) next
  df_ds   <- df_ds[order(df_ds$repetition), ]
  saveRDS(df_ds, file.path(DIR_RESULTS,
          paste0("risultati_ADAC_SWORDSCALE_", safe_ds, "_", start_rep, "_", end_rep, ".rds")))
}
cat("Saved: risultati_finali_benchmark_ADAC/risultati_ADAC_SWORDSCALE_<dataset>_1_50.rds\n")


############################################################
# CYCLE 2 — SWORD perturbed
# foreach over (dataset x rep); saves pred_perturbed_SWORDSCALE_*.rds
############################################################
combos_pert <- expand.grid(
  ds  = names_248datasets,
  rep = seq(start_rep, end_rep_pert),
  stringsAsFactors = FALSE
)
cat("\nCYCLE 2: SWORD perturbed --", nrow(combos_pert), "combinations (ds x rep)\n")

foreach(
  ds_name = combos_pert$ds,
  rep     = combos_pert$rep,
  .options.multicore = list(preschedule = FALSE),
  .errorhandling = "remove"
) %dopar% {
  tryCatch({

    safe_ds  <- gsub("[^A-Za-z0-9_]", "_", ds_name)
    out_file <- file.path(DIR_PERTURB,
                          paste0("pred_perturbed_SWORDSCALE_", safe_ds,
                                 "_rep", rep, ".rds"))
    if (file.exists(out_file)) return(NULL)

    cat("Perturbed:", ds_name, "rep", rep, "\n")

    df <- benchmark_248_dataset_list[[ds_name]]
    colnames(df) <- gsub("-", "_", colnames(df))
    colnames(df) <- gsub(" ", "_", colnames(df))
    if (any(sapply(df, is.factor) | sapply(df, is.character))) {
      df <- fastDummies::dummy_cols(df, remove_selected_columns = TRUE,
                                    remove_first_dummy = TRUE)
    }
    colnames(df) <- sanitize_colnames(colnames(df))
    df <- rename_except_target(df)

    y <- df$target
    x <- df[, setdiff(names(df), "target")]

    if (nrow(df) > 5000) {
      partition_matrix <- matrix(NA, nrow = 3500, ncol = 200)
      set.seed(10)
      for (ii in 1:200) partition_matrix[, ii] <- sample(nrow(df), 3500)
    } else {
      set.seed(10)
      partition_matrix <- caret::createDataPartition(df$target, p = 0.7,
                                                     times = 200, list = FALSE)
    }

    idx_train <- partition_matrix[, rep]
    x_train   <- x[idx_train,  , drop = FALSE]
    y_train   <- y[idx_train]
    x_test    <- x[-idx_train, , drop = FALSE]

    constant_cols <- sapply(x_train, function(col) length(unique(col)) == 1)
    x_train <- x_train[, !constant_cols, drop = FALSE]
    x_test  <- x_test[,  !constant_cols, drop = FALSE]

    # -- Perturb x_train and retrain --
    x_train_pert <- perturb_variables(x_train)

    model_pert <- SVM.ROT.RF.OOB(
      Covariates    = x_train_pert,
      y             = y_train,
      nmin          = 5,
      minleaf       = 2,
      cp            = 0,
      n_perc        = 1,
      n_topCor      = 2,
      threshold_COR = 1,
      cost_C        = 1,
      m             = 100,
      rf_var        = ncol(x_train_pert),
      Weight_Scheme = "scale",
      rand_ntopcor  = TRUE,
      parallel      = FALSE,
      OOB           = FALSE,
      verbose       = FALSE
    )

    # Predict on the ORIGINAL (non-perturbed) x_test
    pred_pert <- rowMeans(SVM.ROT.PRED.RF_all(x_test, model_pert), na.rm = TRUE)
    saveRDS(as.numeric(pred_pert), out_file)
    NULL

  }, error = function(e) {
    cat("Error (pert):", ds_name, "rep", rep, ":", e$message, "\n")
    NULL
  })
}

if (.Platform$OS.type == "windows") stopCluster(cl)
cat("\nEND run_benchmark_SWORD.R\n")
