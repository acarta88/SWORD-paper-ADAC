# ==============================================================================
# FULL BENCHMARK GENERATION — Rotation Forest (ROT) baseline on the 248 real
# datasets. Standard cycle (50 repetitions) + perturbed-variables robustness
# cycle.
#
# NOTE: "from scratch" generation script (slow), included for transparency/
# provenance.
# NOT required to reproduce any paper figure/table.
#
# Rotation Forest has no maintained CRAN package for this use case, so the
# implementation below (RotationForest/predict.RotationForest/BuildOneModel)
# is a from-scratch implementation of Rodriguez et al. 2006, built on rpart.
# k = number of feature groups (rotations); the paper uses
# k = ceiling((n_features)/3), falling back to /5 if PCA on a too-small
# group fails.
# ==============================================================================

library(fastDummies)
library(caret)
library(Metrics)
library(doParallel)
library(foreach)
library(rpart)

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
source(file.path(SCRIPT_DIR, "rotationforest_patch.R"))
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

# fit_rotation_forest() here takes ncol(df) (includes the target column),
# matching this script's original (ncol_df - 1) convention.
fit_rotation_forest_local <- function(x_train, y_train, ncol_df) {
  fit_rotation_forest(x_train, y_train, ncol_df - 1)
}

n_workers <- min(10, parallel::detectCores() - 1L)
if (.Platform$OS.type == "windows") {
  cl <- makeCluster(n_workers)
  registerDoParallel(cl)
  clusterExport(cl, varlist = c("SCRIPT_DIR", "RotationForest", "predict.RotationForest",
                                "BuildOneModel", "fit_rotation_forest", "fit_rotation_forest_local"))
  clusterEvalQ(cl, {
    library(Metrics); library(fastDummies); library(caret); library(rpart)
    source(file.path(SCRIPT_DIR, "benchmark_common.R"))
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
    result_file <- file.path(DIR_RESULTS, paste0("risultati_ROT_ADAC_", safe_ds, "_", start_rep, "_", end_rep, ".rds"))
    if (file.exists(result_file)) { cat("Skip (already done):", ds_name, "\n"); return(NULL) }

    cat("\nStart dataset:", ds_name, "\n")

    df <- preprocess_benchmark_df(benchmark_248_dataset_list[[ds_name]])
    y <- df$target
    x <- df[, setdiff(names(df), "target")]
    partition_matrix <- make_partition_matrix(df)

    risultati_finali <- list()
    for (rep in start_rep:end_rep) {
      idx_train <- partition_matrix[, rep]
      x_train <- x[idx_train, , drop = FALSE]; y_train <- y[idx_train]
      x_test  <- x[-idx_train, , drop = FALSE]; y_test  <- y[-idx_train]

      start_time <- Sys.time()
      set.seed(rep)
      model <- fit_rotation_forest_local(x_train, y_train, ncol(df))
      time_forest <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

      pred_test <- predict.RotationForest(model, x_test, prob = TRUE)

      risultati_finali[[rep]] <- data.frame(
        dataset = ds_name, repetition = rep,
        sd_y_train = sd(y_train, na.rm = TRUE), sd_y_test = sd(y_test, na.rm = TRUE),
        ncolumn = ncol(x),
        RMSE = Metrics::rmse(y_test, pred_test), MAE = Metrics::mae(y_test, pred_test),
        R2 = (cor(y_test, pred_test))^2, RGA = RGA_metric(y_test, pred_test),
        time_forest = time_forest, ntree = 100, mtry = "ROT"
      )

      saveRDS(as.numeric(pred_test), file.path(DIR_PRED, paste0("pred_ROT_", safe_ds, "_rep", rep, ".rds")))
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
      file.exists(file.path(DIR_PERTURB, paste0("pred_perturbed_ROT_", safe_ds, "_rep", r, ".rds")))))
    if (all_done) return(NULL)

    cat("\nPerturbed dataset:", ds_name, "\n")

    df <- preprocess_benchmark_df(benchmark_248_dataset_list[[ds_name]])
    y <- df$target
    x <- df[, setdiff(names(df), "target")]
    partition_matrix <- make_partition_matrix(df)

    for (rep in start_rep:end_rep) {
      out_file <- file.path(DIR_PERTURB, paste0("pred_perturbed_ROT_", safe_ds, "_rep", rep, ".rds"))
      if (file.exists(out_file)) next

      idx_train <- partition_matrix[, rep]
      x_train   <- x[idx_train,  , drop = FALSE]; y_train <- y[idx_train]
      x_test    <- x[-idx_train, , drop = FALSE]

      constant_cols <- sapply(x_train, function(col) length(unique(col)) == 1)
      x_train <- x_train[, !constant_cols, drop = FALSE]
      x_test  <- x_test[,  !constant_cols, drop = FALSE]

      x_train_pert <- perturb_variables(x_train)
      set.seed(rep)
      model_pert <- fit_rotation_forest_local(x_train_pert, y_train, ncol(df))
      pred_pert <- as.numeric(predict.RotationForest(model_pert, x_test, prob = TRUE))
      saveRDS(pred_pert, out_file)
    }
    gc(); NULL
  }, error = function(e) {
    cat("Error (perturbed):", ds_name, "-", e$message, "\n")
    NULL
  })
}

if (.Platform$OS.type == "windows") stopCluster(cl)
cat("\nEND run_benchmark_ROT.R\n")
