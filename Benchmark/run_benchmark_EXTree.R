# ==============================================================================
# FULL BENCHMARK GENERATION — Extremely Randomized Trees (EXTree) baseline on
# the 248 real datasets. Standard cycle (50 repetitions) + perturbed-variables
# robustness cycle.
#
# NOTE: "from scratch" generation script (slow), included for transparency/
# provenance.
# NOT required to reproduce any paper figure/table.
#
# Recipe: extraTrees::extraTrees(x_train, y_train, ntree = 100).
# WARNING: the `extraTrees` package is a thin R wrapper around a Java
# implementation and requires `rJava` + a working JDK on the machine. It is
# no longer actively maintained on CRAN; if installation fails, this is a
# dependency problem with `extraTrees`/`rJava`, not with this script.
# ==============================================================================

library(fastDummies)
library(caret)
library(Metrics)
options(java.parameters = "-Xmx8g")
library(rJava)
library(extraTrees)

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

############################################################
# CYCLE 1 — standard
############################################################
for (ds_name in names_248datasets) {
  safe_ds <- gsub("[^A-Za-z0-9_]", "_", ds_name)
  result_file <- file.path(DIR_RESULTS, paste0("risultati_EXTree_ADAC_", safe_ds, "_", start_rep, "_", end_rep, ".rds"))
  if (file.exists(result_file)) { cat("Skip (already done):", ds_name, "\n"); next }

  tryCatch({
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

      constant_cols <- sapply(x_train, function(col) length(unique(col)) == 1)
      x_train <- x_train[, !constant_cols, drop = FALSE]
      x_test  <- x_test[,  !constant_cols, drop = FALSE]

      start_time <- Sys.time()
      set.seed(rep)
      model <- extraTrees(x_train, y_train, ntree = 100)
      time_forest <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

      pred_matrix_all <- predict(model, x_test, allValues = TRUE)
      pred_test <- rowMeans(pred_matrix_all, na.rm = TRUE)

      cor_matrix    <- cor(pred_matrix_all, use = "pairwise.complete.obs")
      cor_tree      <- mean(cor_matrix[upper.tri(cor_matrix)], na.rm = TRUE)
      r2_trees      <- apply(pred_matrix_all, 2, function(p) (cor(y_test, p, use = "pairwise.complete.obs"))^2)
      strength_tree <- mean(r2_trees, na.rm = TRUE)

      saveRDS(as.numeric(pred_test), file.path(DIR_PRED, paste0("pred_EXTree_", safe_ds, "_rep", rep, ".rds")))

      risultati_finali[[rep]] <- data.frame(
        dataset = ds_name, repetition = rep,
        sd_y_train = sd(y[idx_train], na.rm = TRUE), sd_y_test = sd(y_test, na.rm = TRUE),
        ncolumn = ncol(x),
        RMSE = Metrics::rmse(y_test, pred_test), MAE = Metrics::mae(y_test, pred_test),
        R2 = (cor(y_test, pred_test))^2, RGA = RGA_metric(y_test, pred_test),
        cor_tree = cor_tree, strength_tree = strength_tree,
        time_forest = time_forest, mtry = "EXTree", ntree = 100
      )
    }

    saveRDS(do.call(rbind, risultati_finali), file = result_file)
    gc()
  }, error = function(e) {
    cat("Error processing dataset:", ds_name, "-", e$message, "\n")
  })
}

############################################################
# CYCLE 2 — perturbed
############################################################
for (ds_name in names_248datasets) {
  safe_ds <- gsub("[^A-Za-z0-9_]", "_", ds_name)
  all_done <- all(sapply(start_rep:end_rep, function(r)
    file.exists(file.path(DIR_PERTURB, paste0("pred_perturbed_EXTree_", safe_ds, "_rep", r, ".rds")))))
  if (all_done) next

  tryCatch({
    cat("\nPerturbed dataset:", ds_name, "\n")

    df <- preprocess_benchmark_df(benchmark_248_dataset_list[[ds_name]])
    y <- df$target
    x <- df[, setdiff(names(df), "target")]
    partition_matrix <- make_partition_matrix(df)

    for (rep in start_rep:end_rep) {
      out_file <- file.path(DIR_PERTURB, paste0("pred_perturbed_EXTree_", safe_ds, "_rep", rep, ".rds"))
      if (file.exists(out_file)) next

      idx_train <- partition_matrix[, rep]
      x_train   <- x[idx_train,  , drop = FALSE]; y_train <- y[idx_train]
      x_test    <- x[-idx_train, , drop = FALSE]

      constant_cols <- sapply(x_train, function(col) length(unique(col)) == 1)
      x_train <- x_train[, !constant_cols, drop = FALSE]
      x_test  <- x_test[,  !constant_cols, drop = FALSE]

      x_train_pert <- perturb_variables(x_train)
      set.seed(rep)
      model_pert <- extraTrees(x_train_pert, y_train, ntree = 100)
      pred_pert  <- rowMeans(predict(model_pert, x_test, allValues = TRUE), na.rm = TRUE)
      saveRDS(as.numeric(pred_pert), out_file)
    }
    gc()
  }, error = function(e) {
    cat("Error (perturbed):", ds_name, "-", e$message, "\n")
  })
}

cat("\nEND run_benchmark_EXTree.R\n")
