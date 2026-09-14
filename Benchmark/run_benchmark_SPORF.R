# ==============================================================================
# FULL BENCHMARK GENERATION — Sparse Projection Oblique Randomer Forest
# (SPORF) baseline on the 248 real datasets. Standard cycle (50 repetitions)
# + perturbed-variables robustness cycle.
#
# NOTE: "from scratch" generation script (slow), included for transparency/
# provenance.
# NOT required to reproduce any paper figure/table.
#
# Recipe: same ODRF::ODRF() call as run_benchmark_ODRF.R, but with
# NodeRotateFun = "RotMatRand" (random sparse rotations), which approximates
# the SPORF methodology within the ODRF package. Uses the same hand-written
# predict methods (odrf_predict_patch.R). NOTE: unlike every other model in
# this benchmark, the original source for SPORF does not call
# rename_except_target() on the predictors -- preserved here for fidelity.
# ==============================================================================

library(fastDummies)
library(caret)
library(Metrics)
library(doParallel)
library(foreach)
library(ODRF)

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

source(file.path(SCRIPT_DIR, "benchmark_common.R"))
source(file.path(SCRIPT_DIR, "odrf_predict_patch.R"))
load(file.path(DATA, "benchmark_248_dataset_list_ADAC.rda"))
names_248datasets <- names(benchmark_248_dataset_list)

DIR_RESULTS <- file.path(SCRIPT_DIR, "risultati_finali_benchmark_ADAC")
DIR_PRED    <- file.path(SCRIPT_DIR, "vettori_predizioni_benchmark_ADAC")
DIR_PERTURB <- file.path(SCRIPT_DIR, "vettori_perturbed_benchmark_ADAC")
dir.create(DIR_RESULTS,  showWarnings = FALSE)
dir.create(DIR_PRED,     showWarnings = FALSE)
dir.create(DIR_PERTURB,  showWarnings = FALSE)

start_rep <- 1
end_rep   <- 50

n_workers <- min(10, parallel::detectCores() - 1L)
if (.Platform$OS.type == "windows") {
  cl <- makeCluster(n_workers)
  registerDoParallel(cl)
  clusterExport(cl, varlist = c("PROJ", "SCRIPT_DIR"))
  clusterEvalQ(cl, {
    library(Metrics); library(fastDummies); library(caret); library(ODRF)
    source(file.path(SCRIPT_DIR, "benchmark_common.R"))
    source(file.path(SCRIPT_DIR, "odrf_predict_patch.R"))
  })
  clusterExport(cl, varlist = c("benchmark_248_dataset_list", "DIR_RESULTS", "DIR_PRED"))
} else {
  library(doMC)
  registerDoMC(cores = n_workers)
}

############################################################
# CYCLE 1 — standard
############################################################
foreach(
  ds_name = names_248datasets,
  .errorhandling = "remove",
  .options.multicore = list(preschedule = FALSE)
) %dopar% {
  tryCatch({
    safe_ds <- gsub("[^A-Za-z0-9_]", "_", ds_name)
    result_file <- file.path(DIR_RESULTS, paste0("risultati_SPORF_ADAC_", safe_ds, "_", start_rep, "_", end_rep, ".rds"))
    if (file.exists(result_file)) { cat("Skip (already done):", ds_name, "\n"); return(NULL) }

    cat("\nStart dataset:", ds_name, "\n")

    # NOTE: no rename_except_target() here -- see header note.
    df <- preprocess_benchmark_df(benchmark_248_dataset_list[[ds_name]], do_rename = FALSE)
    y <- df$target
    x <- df[, setdiff(names(df), "target")]
    partition_matrix <- make_partition_matrix(df)

    risultati_finali <- list()
    for (rep in start_rep:end_rep) {
      idx_train <- partition_matrix[, rep]
      x_train <- x[idx_train, , drop = FALSE]; y_train <- y[idx_train]
      x_test  <- x[-idx_train, , drop = FALSE]; y_test  <- y[-idx_train]

      df_train <- data.frame(x_train, target = y_train)
      df_test  <- data.frame(x_test,  target = y_test)
      constant_cols <- sapply(df_train, function(col) length(unique(col)) == 1)
      df_train <- df_train[, !constant_cols]
      df_test  <- df_test[,  !constant_cols]

      start_time <- Sys.time()
      set.seed(10)
      model <- ODRF(target ~ ., data = df_train, split = "mse", parallel = FALSE,
                    ntrees = 100, storeOOB = FALSE, NodeRotateFun = "RotMatRand")
      time_forest <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

      pred_matrix_all <- predict.ODRF(model, df_test[, -which(names(df_test) == "target")], type = "tree")
      pred_test <- rowMeans(pred_matrix_all, na.rm = TRUE)

      cor_matrix    <- cor(pred_matrix_all, use = "pairwise.complete.obs")
      cor_tree      <- mean(cor_matrix[upper.tri(cor_matrix)], na.rm = TRUE)
      r2_trees      <- apply(pred_matrix_all, 2, function(p) (cor(y_test, p, use = "pairwise.complete.obs"))^2)
      strength_tree <- mean(r2_trees, na.rm = TRUE)

      saveRDS(as.numeric(pred_test), file.path(DIR_PRED, paste0("pred_SPORF_", safe_ds, "_rep", rep, ".rds")))

      risultati_finali[[rep]] <- data.frame(
        dataset = ds_name, repetition = rep,
        sd_y_train = sd(y[idx_train], na.rm = TRUE), sd_y_test = sd(y_test, na.rm = TRUE),
        ncolumn = ncol(x),
        RMSE = Metrics::rmse(y_test, pred_test), MAE = Metrics::mae(y_test, pred_test),
        R2 = (cor(y_test, pred_test))^2, RGA = RGA_metric(y_test, pred_test),
        cor_tree = cor_tree, strength_tree = strength_tree,
        time_forest = time_forest, mtry = "SPORF", ntree = 100
      )
    }

    saveRDS(do.call(rbind, risultati_finali), file = result_file)
    gc(); NULL
  }, error = function(e) {
    cat("Error processing dataset:", ds_name, "-", e$message, "\n")
    NULL
  })
}

############################################################
# CYCLE 2 — perturbed
############################################################
foreach(
  ds_name = names_248datasets,
  .errorhandling = "remove",
  .options.multicore = list(preschedule = FALSE)
) %dopar% {
  tryCatch({
    safe_ds <- gsub("[^A-Za-z0-9_]", "_", ds_name)
    all_done <- all(sapply(start_rep:end_rep, function(r)
      file.exists(file.path(DIR_PERTURB, paste0("pred_perturbed_SPORF_", safe_ds, "_rep", r, ".rds")))))
    if (all_done) return(NULL)

    cat("\nPerturbed dataset:", ds_name, "\n")

    df <- preprocess_benchmark_df(benchmark_248_dataset_list[[ds_name]], do_rename = FALSE)
    y <- df$target
    x <- df[, setdiff(names(df), "target")]
    partition_matrix <- make_partition_matrix(df)

    for (rep in start_rep:end_rep) {
      out_file <- file.path(DIR_PERTURB, paste0("pred_perturbed_SPORF_", safe_ds, "_rep", rep, ".rds"))
      if (file.exists(out_file)) next

      idx_train <- partition_matrix[, rep]
      x_train   <- x[idx_train,  , drop = FALSE]; y_train <- y[idx_train]
      x_test    <- x[-idx_train, , drop = FALSE]

      df_train <- data.frame(x_train, target = y_train)
      df_test  <- data.frame(x_test,  target = rep(NA, nrow(x_test)))
      constant_cols <- sapply(df_train, function(col) length(unique(col)) == 1)
      df_train <- df_train[, !constant_cols]
      df_test  <- df_test[,  !constant_cols]

      x_train_clean <- df_train[, setdiff(names(df_train), "target")]
      x_test_clean  <- df_test[,  setdiff(names(df_test),  "target")]
      x_train_pert  <- perturb_variables(x_train_clean)
      df_train_pert <- data.frame(x_train_pert, target = y_train)

      set.seed(rep)
      model_pert <- ODRF(target ~ ., data = df_train_pert, split = "mse", parallel = FALSE,
                         ntrees = 100, storeOOB = FALSE, NodeRotateFun = "RotMatRand")
      pred_pert <- rowMeans(predict.ODRF(model_pert, x_test_clean, type = "tree"), na.rm = TRUE)
      saveRDS(as.numeric(pred_pert), out_file)
    }
    gc(); NULL
  }, error = function(e) {
    cat("Error (perturbed):", ds_name, "-", e$message, "\n")
    NULL
  })
}

if (.Platform$OS.type == "windows") stopCluster(cl)
cat("\nEND run_benchmark_SPORF.R\n")
