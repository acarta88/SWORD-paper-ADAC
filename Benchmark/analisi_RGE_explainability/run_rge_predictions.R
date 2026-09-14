# ==============================================================================
# FULL RGE GENERATION — leave-one-variable-out predictions for all 11 models
# (SWORD + 10 baselines) on the 248 real datasets.
#
# For every (dataset, variable) pair with >= 3 predictors, the variable is
# removed, every model is retrained on the remaining predictors, and its
# test-set predictions are saved. These per-task prediction vectors are later
# compared against the "full-model" predictions (from the main benchmark) to
# compute the RGE (Rank Graduation Explainability) score used in
# fig:rge_vs_sword -- a rank-based measure of how much removing a variable
# changes what the model predicts, i.e. how much the model actually "uses"
# that variable.
#
# NOTE: "from scratch" generation script (VERY slow: 248 datasets x ~10-20
# variables each x 11 models). Included for full transparency/provenance of
# the results already shipped in pred_RGE_benchmark/ (raw, ~8.3 GB, not
# deposited in this repository) and output_RGE/df_RGE_wide.rds
# (aggregated, 0.3 MB, deposited). NOT required to reproduce fig:rge_vs_sword
# -- fig6_rge_explainability.R reads directly from output_RGE/df_RGE_wide.rds.
#
# Only 1 repetition per task (start_rep = end_rep = 1): the train/test split
# used is repetition 1 of the same 200-partition matrix as the main
# benchmark (same set.seed(10) + createDataPartition/sample recipe).
# ==============================================================================

library(fastDummies)
library(caret)
library(Metrics)
library(doParallel)
library(foreach)
library(randomForest)
library(RRF)
library(aorsf)
library(partykit)
library(rpart)
library(ODRF)
library(e1071)
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
SCRIPT_DIR <- .get_script_dir()          # .../Benchmark/analisi_RGE_explainability
BENCH_DIR  <- normalizePath(file.path(SCRIPT_DIR, ".."))   # .../Benchmark
PROJ       <- normalizePath(file.path(BENCH_DIR, ".."))    # repo root
DATA       <- file.path(PROJ, "data")

source(file.path(PROJ, "R", "SVM.ROT_Functions_ADAC.R"))
source(file.path(BENCH_DIR, "benchmark_common.R"))
source(file.path(BENCH_DIR, "odrf_predict_patch.R"))
source(file.path(BENCH_DIR, "rotationforest_patch.R"))

load(file.path(DATA, "benchmark_248_dataset_list_ADAC.rda"))
names_248datasets <- names(benchmark_248_dataset_list)

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

out_dir <- file.path(SCRIPT_DIR, "pred_RGE_benchmark")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

start_rep <- 1
end_rep   <- 1

########## Data-prep helpers specific to the leave-one-variable-out design ##########

# Order for each (dataset, variable) task: sanitize raw column names -> drop
# the target variable -> dummy-encode -> rename predictors to varA/varB/...
prepare_data <- function(ds_name, var_orig, idx_train) {
  df_raw <- benchmark_248_dataset_list[[ds_name]]
  colnames(df_raw) <- sanitize_colnames(colnames(df_raw))

  df <- df_raw[, setdiff(colnames(df_raw), var_orig), drop = FALSE]
  if (any(sapply(df, is.factor) | sapply(df, is.character)))
    df <- fastDummies::dummy_cols(df, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
  df <- rename_except_target(df)

  y       <- df$target
  x       <- df[, setdiff(names(df), "target"), drop = FALSE]
  x_train <- x[ idx_train, , drop = FALSE]
  y_train <- y[ idx_train]
  x_test  <- x[-idx_train, , drop = FALSE]

  constant_cols <- sapply(x_train, function(col) length(unique(col)) == 1)
  x_train <- x_train[, !constant_cols, drop = FALSE]
  x_test  <- x_test[,  !constant_cols, drop = FALSE]

  list(x_train = x_train, y_train = y_train, x_test = x_test)
}

# Same, but returns a formula-ready df_train (for ODRF/aorsf/CF/SPORF).
prepare_data_df <- function(ds_name, var_orig, idx_train) {
  df_raw <- benchmark_248_dataset_list[[ds_name]]
  colnames(df_raw) <- sanitize_colnames(colnames(df_raw))

  df <- df_raw[, setdiff(colnames(df_raw), var_orig), drop = FALSE]
  if (any(sapply(df, is.factor) | sapply(df, is.character)))
    df <- fastDummies::dummy_cols(df, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
  df <- rename_except_target(df)

  y        <- df$target
  x        <- df[, setdiff(names(df), "target"), drop = FALSE]
  df_train <- data.frame(x[idx_train, , drop = FALSE], target = y[idx_train])
  x_test   <- x[-idx_train, , drop = FALSE]

  constant_cols <- sapply(df_train, function(col) length(unique(col)) == 1)
  df_train <- df_train[, !constant_cols, drop = FALSE]
  x_test   <- x_test[, intersect(setdiff(colnames(df_train), "target"), colnames(x_test)), drop = FALSE]

  list(df_train = df_train, x_test = x_test)
}

# set.seed(10) is fixed -> identical partition matrix for the same dataset
# regardless of which variable is being removed or which worker computes it.
get_partition <- function(ds_name) {
  df_raw <- benchmark_248_dataset_list[[ds_name]]
  colnames(df_raw) <- sanitize_colnames(colnames(df_raw))
  df_full <- df_raw
  if (any(sapply(df_full, is.factor) | sapply(df_full, is.character)))
    df_full <- fastDummies::dummy_cols(df_full, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
  df_full <- rename_except_target(df_full)
  make_partition_matrix(df_full)
}

########## Task list (dataset x variable), sorted by increasing #variables ##########

tasks <- do.call(rbind, lapply(names_248datasets, function(ds_name) {
  df_raw <- benchmark_248_dataset_list[[ds_name]]
  colnames(df_raw) <- sanitize_colnames(colnames(df_raw))
  vars <- setdiff(colnames(df_raw), "target")
  if (length(vars) < 3) return(NULL)  # removing any var would leave < 2 predictors (pre-dummy)
  data.frame(ds_name = ds_name, var_orig = vars, stringsAsFactors = FALSE)
}))
nvars_per_ds    <- table(tasks$ds_name)[tasks$ds_name]
tasks           <- tasks[order(as.integer(nvars_per_ds)), ]
rownames(tasks) <- NULL
task_list       <- split(tasks, seq_len(nrow(tasks)))

########## Parallel backend ##########

n_workers <- min(10, parallel::detectCores() - 1L)
if (.Platform$OS.type == "windows") {
  cl <- makeCluster(n_workers)
  registerDoParallel(cl)
  clusterExport(cl, varlist = c("PROJ", "BENCH_DIR", "out_dir", "start_rep", "end_rep"))
  clusterEvalQ(cl, {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    library(fastDummies); library(caret);    library(Metrics)
    library(randomForest); library(RRF);    library(aorsf)
    library(partykit);    library(rpart);   library(ODRF)
    library(e1071);        library(extraTrees)
    source(file.path(PROJ, "R", "SVM.ROT_Functions_ADAC.R"))
    source(file.path(BENCH_DIR, "benchmark_common.R"))
    source(file.path(BENCH_DIR, "odrf_predict_patch.R"))
    source(file.path(BENCH_DIR, "rotationforest_patch.R"))
  })
  clusterExport(cl, varlist = c("benchmark_248_dataset_list", "prepare_data",
                                "prepare_data_df", "get_partition"))
} else {
  library(doMC)
  registerDoMC(cores = n_workers)
}

########## Helper: run one (task x model) leave-one-var-out prediction ##########

run_rge_task <- function(model_name, fit_predict_fn) {
  foreach(
    task = task_list,
    .options.multicore = list(preschedule = FALSE),
    .errorhandling = "remove"
  ) %dopar% {
    tryCatch({
      ds_name  <- task$ds_name
      var_orig <- task$var_orig
      safe_ds  <- gsub("[^A-Za-z0-9_]", "_", ds_name)
      safe_var <- substr(gsub("[^A-Za-z0-9_]", "_", var_orig), 1, 40)

      pm <- get_partition(ds_name)

      for (rep in start_rep:end_rep) {
        out_file <- file.path(out_dir, paste0("pred_RGE_", model_name, "_", safe_ds,
                                              "_VARIAB_", safe_var, "_rep", rep, ".rds"))
        if (file.exists(out_file)) next
        pred <- fit_predict_fn(ds_name, var_orig, pm[, rep], rep)
        saveRDS(as.numeric(pred), out_file)
      }
    }, error = function(e) cat("Error", model_name, ":", task$ds_name, task$var_orig, "-", e$message, "\n"))
  }
}

############################################################
# SWORD
############################################################
run_rge_task("SWORD", function(ds_name, var_orig, idx_train, rep) {
  d <- prepare_data(ds_name, var_orig, idx_train)
  set.seed(rep)
  model <- SVM.ROT.RF.OOB(
    Covariates = d$x_train, y = d$y_train,
    nmin = 5, minleaf = 2, cp = 0, n_perc = 1,
    n_topCor = 2, threshold_COR = 1, Weight_Scheme = "scale", relation = "Pearson",
    cost_C = 1, m = 100, rf_var = ncol(d$x_train),
    rand_ntopcor = TRUE, parallel = FALSE, OOB = FALSE, verbose = FALSE
  )
  rowMeans(SVM.ROT.PRED.RF_all(d$x_test, model), na.rm = TRUE)
})

############################################################
# Random Forest
############################################################
run_rge_task("RF", function(ds_name, var_orig, idx_train, rep) {
  d <- prepare_data(ds_name, var_orig, idx_train)
  set.seed(rep)
  model <- randomForest::randomForest(d$x_train, d$y_train, ntree = 100)
  predict(model, d$x_test)
})

############################################################
# ODRF
############################################################
run_rge_task("ODRF", function(ds_name, var_orig, idx_train, rep) {
  d <- prepare_data_df(ds_name, var_orig, idx_train)
  set.seed(rep)
  model <- ODRF::ODRF(target ~ ., data = d$df_train, split = "mse", parallel = FALSE,
                      ntrees = 100, storeOOB = FALSE)
  rowMeans(predict.ODRF(model, d$x_test, type = "tree"), na.rm = TRUE)
})

############################################################
# Rotation Forest (ROT)
############################################################
run_rge_task("ROT", function(ds_name, var_orig, idx_train, rep) {
  d <- prepare_data(ds_name, var_orig, idx_train)
  set.seed(rep)
  model <- fit_rotation_forest(d$x_train, d$y_train, ncol(d$x_train))
  predict.RotationForest(model, d$x_test)
})

############################################################
# Extra Trees (EXTree)
############################################################
run_rge_task("EXTree", function(ds_name, var_orig, idx_train, rep) {
  d <- prepare_data(ds_name, var_orig, idx_train)
  set.seed(rep)
  model <- extraTrees::extraTrees(d$x_train, d$y_train, ntree = 100)
  rowMeans(predict(model, d$x_test, allValues = TRUE), na.rm = TRUE)
})

############################################################
# RRF (Regularized Random Forest)
############################################################
run_rge_task("RRF", function(ds_name, var_orig, idx_train, rep) {
  d <- prepare_data(ds_name, var_orig, idx_train)
  set.seed(rep)
  model <- RRF::RRF(d$x_train, d$y_train, ntree = 100)
  rowMeans(predict(model, d$x_test, predict.all = TRUE)$individual, na.rm = TRUE)
})

############################################################
# aorsf (GLM node models)
############################################################
run_rge_task("aorsf", function(ds_name, var_orig, idx_train, rep) {
  d <- prepare_data_df(ds_name, var_orig, idx_train)
  set.seed(rep)
  model <- aorsf::orsf(d$df_train, target ~ ., n_tree = 100,
                       control = aorsf::orsf_control_regression(method = "glm"))
  rowMeans(predict(model, d$x_test, pred_type = "mean", pred_aggregate = FALSE), na.rm = TRUE)
})

############################################################
# aorsfNET (penalized node models)
############################################################
run_rge_task("aorsfNET", function(ds_name, var_orig, idx_train, rep) {
  d <- prepare_data_df(ds_name, var_orig, idx_train)
  set.seed(rep)
  model <- aorsf::orsf(d$df_train, target ~ ., n_tree = 100,
                       control = aorsf::orsf_control_regression(method = "net"))
  rowMeans(predict(model, d$x_test, pred_type = "mean", pred_aggregate = FALSE), na.rm = TRUE)
})

############################################################
# SPORF (ODRF with NodeRotateFun = "RotMatRand")
############################################################
run_rge_task("SPORF", function(ds_name, var_orig, idx_train, rep) {
  d <- prepare_data_df(ds_name, var_orig, idx_train)
  set.seed(10)
  model <- ODRF::ODRF(target ~ ., data = d$df_train, split = "mse", NodeRotateFun = "RotMatRand",
                      parallel = FALSE, ntrees = 100, storeOOB = FALSE)
  rowMeans(predict.ODRF(model, d$x_test, type = "tree"), na.rm = TRUE)
})

############################################################
# Conditional Forest (CF)
############################################################
run_rge_task("CF", function(ds_name, var_orig, idx_train, rep) {
  d <- prepare_data_df(ds_name, var_orig, idx_train)
  set.seed(10)
  model      <- partykit::cforest(target ~ ., data = d$df_train, ntree = 100)
  batch_size <- 2000
  n_rows     <- nrow(d$x_test)
  pred_test  <- numeric(n_rows)
  for (i in seq_len(ceiling(n_rows / batch_size))) {
    si <- (i - 1) * batch_size + 1; ei <- min(i * batch_size, n_rows)
    pred_test[si:ei] <- predict(model, d$x_test[si:ei, , drop = FALSE])
  }
  pred_test
})

############################################################
# Linear SVM (SVMlin)
############################################################
run_rge_task("SVMlin", function(ds_name, var_orig, idx_train, rep) {
  d <- prepare_data(ds_name, var_orig, idx_train)
  set.seed(rep)
  model <- e1071::svm(d$x_train, d$y_train, kernel = "linear", type = "eps-regression")
  predict(model, d$x_test)
})

if (.Platform$OS.type == "windows") {
  stopCluster(cl)
  registerDoSEQ()
}

cat("\nEND run_rge_predictions.R\n")
