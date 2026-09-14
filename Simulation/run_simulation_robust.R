# ==============================================================================
# FULL SIMULATION — Weight_Scheme = "robust", all seeds
# Parallelized over DATASETS (foreach/doMC on Linux, makeCluster on Windows)
# Serial 36-combo grid inside each worker
# Output: risultati_full_simulation_robust/
#
# NOTE: this is the "from scratch" generation script (slow: hours, meant for a
# many-core machine), included for full transparency/provenance. Its output
# folder (risultati_full_simulation_robust/) is not shipped in this repo; the
# light cache that IS shipped (results_cache/df_SWORD72_cached.rds) already
# includes the robust-scheme results the paper's figures/tables need (see
# tab1_wilcox_nrmse.R, tab3_tabA1_wilcox_corstrength.R, fig2_tabA2_compute_time.R,
# fig1_variable_importance.R) — this script is NOT required to reproduce any of them. Run
# this only to regenerate everything from the raw simulated datasets.
# ==============================================================================
rm(list = ls()); gc()

# ==============================================================================
# PART 0: SETUP
# ==============================================================================

pkgs_needed <- c(
  "data.table", "infotheo", "WeightSVM", "data.tree",
  "stringr", "Metrics", "fastDummies",
  "foreach", "doParallel"
)
pkgs_missing <- pkgs_needed[!sapply(pkgs_needed, requireNamespace, quietly = TRUE)]
if (length(pkgs_missing) > 0) {
  cat("Installing:", paste(pkgs_missing, collapse = ", "), "\n")
  install.packages(pkgs_missing, repos = "https://cloud.r-project.org")
}

library(data.table)
library(infotheo)
library(WeightSVM)
library(data.tree)
library(stringr)
library(Metrics)
library(fastDummies)
library(foreach)
library(doParallel)

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

# --- Source core SWORD functions (the only external file needed) ---
source(file.path(PROJ, "R", "SVM.ROT_Functions_ADAC.R"))

# --- Load simulated datasets (all seeds) ---
cat("Loading simulated datasets...\n")
load(file.path(DATA, "simulated_datasets_ADAC.rda"))
# list_simulated_df contains all seeds -- do NOT filter
cat("Available seeds:", paste(names(list_simulated_df), collapse = ", "), "\n")


# ==============================================================================
# PART 1: HELPER FUNCTIONS (inline -- no source of SWORD_simulated_functions)
# ==============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

correlazione_alberi_rf <- function(OOB_matrix) {
  m <- ncol(OOB_matrix)
  if (m < 2) return(NA_real_)
  vals <- numeric(0)
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      ok <- !is.na(OOB_matrix[, i]) & !is.na(OOB_matrix[, j])
      if (sum(ok) > 2) {
        xi <- OOB_matrix[ok, i]; xj <- OOB_matrix[ok, j]
        if (sd(xi) > 0 && sd(xj) > 0)
          vals <- c(vals, cor(xi, xj, use = "pairwise.complete.obs"))
      }
    }
  }
  if (length(vals) == 0) NA_real_ else mean(vals, na.rm = TRUE)
}

strength_rf_breiman <- function(OOB_matrix, target) {
  T <- ncol(OOB_matrix)
  sapply(seq_len(T), function(t) {
    ok <- !is.na(OOB_matrix[, t])
    if (any(ok)) cor(target[ok], OOB_matrix[ok, t],
                     use = "pairwise.complete.obs")^2 else NA_real_
  }) |> mean(na.rm = TRUE)
}

RGA_metric <- function(y_true, y_pred) {
  y_true <- as.numeric(y_true)
  y_pred <- as.numeric(y_pred)
  ok     <- is.finite(y_true) & is.finite(y_pred)
  y_true <- y_true[ok]
  y_pred <- y_pred[ok]
  n <- length(y_true)
  if (n < 2) return(NA_real_)
  I      <- (1:n) / n
  r_pred <- rank(y_pred, ties.method = "min")
  # for ties in y_pred, replace y with the group mean
  y_support <- unname(tapply(y_true, r_pred, mean)[as.character(r_pred)])
  y_sorted  <- y_support[order(r_pred)]
  conc <- 2 * sum(I * y_sorted)
  inc  <- 2 * sum(I * sort(y_true))
  dec  <- 2 * sum(I * sort(y_true, decreasing = TRUE))
  if (inc == dec) return(NA_real_)
  (conc - dec) / (inc - dec)
}

# Flattens the nested list into a list of list(path, df, n_obs, features)
flatten_lista <- function(lst, path = "") {
  if (is.data.frame(lst)) {
    # Extract n_obs and features from the path so we can sort by size
    el <- strsplit(path, "\\$")[[1]]
    el <- el[nzchar(el)]
    n_obs    <- if (length(el) >= 2) as.numeric(gsub("_n_obs", "", el[2])) else 0
    features <- if (length(el) >= 3) as.numeric(gsub("_n_feat", "", el[3])) else 0
    return(list(list(path = path, df = lst, n_obs = n_obs, features = features)))
  }
  if (is.list(lst)) {
    out <- list()
    nomi <- names(lst)
    for (j in seq_along(lst)) {
      nome <- if (!is.null(nomi) && nzchar(nomi[j])) nomi[j] else paste0("[[", j, "]]")
      out <- c(out, flatten_lista(lst[[j]], paste0(path, "$", nome)))
    }
    return(out)
  }
  list()
}

# Parses metadata from the path (7 levels: seed, n_obs, n_feat, Noise, NoLin, cat, err)
parse_path_info <- function(path) {
  el <- strsplit(path, "\\$")[[1]]
  el <- el[nzchar(el)]
  if (length(el) == 7) {
    data.frame(
      seed     = as.numeric(gsub("\\D",    "", el[1])),
      n_obs    = as.numeric(gsub("_n_obs", "", el[2])),
      features = as.numeric(gsub("_n_feat","", el[3])),
      noise    = as.numeric(gsub("_Noise", "", el[4])),
      nonlin   = as.numeric(gsub("_NoLin", "", el[5])),
      cat      = as.numeric(gsub("_cat",   "", el[6])),
      error_sd = as.numeric(gsub("_err",   "", el[7]))
    )
  } else {
    data.frame(seed=NA, n_obs=NA, features=NA, noise=NA, nonlin=NA, cat=NA, error_sd=NA)
  }
}


# ==============================================================================
# PART 2: check_progress() -- call from the console while it is running
# ==============================================================================

check_progress <- function(cartella = file.path(SCRIPT_DIR, "risultati_full_simulation_robust")) {
  rds <- list.files(cartella, pattern = "\\.rds$")
  rds <- rds[rds != "risultati_completi.rds"]
  n_done <- length(rds)

  log_file <- file.path(cartella, "progress_log.txt")
  if (!file.exists(log_file)) {
    cat("Simulation not started yet.\n"); return(invisible(n_done))
  }
  log_lines <- readLines(log_file, warn = FALSE)
  n_totali <- NA
  if (length(log_lines) > 0 && grepl("^TOTALE:", log_lines[1])) {
    n_totali <- as.integer(sub("TOTALE: ", "", log_lines[1]))
    log_lines <- log_lines[-1]
  }

  cat("============================================\n")
  if (!is.na(n_totali)) {
    pct <- round(100 * n_done / n_totali, 1)
    cat(sprintf("  Progress: %d / %d (%.1f%%)\n", n_done, n_totali, pct))
    done_lines <- log_lines[grepl("^\\d{4}-", log_lines)]
    if (length(done_lines) >= 2) {
      t1 <- as.POSIXct(substr(done_lines[1], 1, 19), format = "%Y-%m-%d %H:%M:%S")
      t2 <- as.POSIXct(substr(tail(done_lines,1), 1, 19), format = "%Y-%m-%d %H:%M:%S")
      el <- as.numeric(difftime(t2, t1, units = "mins"))
      if (el > 0 && n_done > 0) {
        rate <- n_done / el
        rem  <- (n_totali - n_done) / rate
        cat(sprintf("  Rate: %.1f datasets/min\n", rate))
        cat(sprintf("  Remaining: %.0f min (%.1f h)\n", rem, rem / 60))
      }
    }
  } else {
    cat("  Datasets completed:", n_done, "\n")
  }
  cat("============================================\n")
  done_lines <- log_lines[grepl("^\\d{4}-", log_lines)]
  if (length(done_lines) > 0) {
    cat("\nLast completed:\n")
    cat(paste(" ", tail(done_lines, 10), collapse = "\n"), "\n")
  }
  invisible(n_done)
}


# ==============================================================================
# PART 3: CONFIGURATION
# ==============================================================================

n_core_totali   <- parallel::detectCores()
# Adaptive worker count based on dataset size:
#   large (n_obs*features >= 5000): few workers to keep RAM in check
#   small: max workers
n_workers_large <- min(70, max(1, n_core_totali - 16))
n_workers_small <- min(70, max(1, n_core_totali - 16))
cartella_robust <- file.path(SCRIPT_DIR, "risultati_full_simulation_robust")
log_file        <- file.path(cartella_robust, "progress_log.txt")
checkpoint_file <- file.path(cartella_robust, "checkpoint.txt")

# rf_var fractions and fixed grid parameters
rf_var_fracs <- c(0.3, 0.6, 1.0)
PARAM_RELATION     <- c("Pearson", "MI")
PARAM_RAND_NTOPCOR <- c(TRUE, FALSE)
PARAM_COST_C       <- c(0.1, 1, 10)

cat("============================================================\n")
cat("  SWORD SIMULATION -- Weight_Scheme = robust\n")
cat("  Platform:", .Platform$OS.type, "\n")
cat("  Cores:", n_core_totali,
    "| Large workers:", n_workers_large,
    "| Small workers:", n_workers_small, "\n")
cat("  Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n\n")

# ==============================================================================
# PART 4: FLATTEN + SORT by size (n_obs DESC, features DESC)
# ==============================================================================

cat("Flattening dataset list...\n")
flat <- flatten_lista(list_simulated_df)
n_totali <- length(flat)
cat("Total datasets (all seeds):", n_totali, "\n")

# Sort: largest first (n_obs DESC, features DESC)
ord <- order(
  sapply(flat, `[[`, "n_obs")    * (-1),
  sapply(flat, `[[`, "features") * (-1)
)
flat <- flat[ord]

cat("Execution order (first 5):\n")
for (i in 1:min(5, length(flat))) {
  cat(sprintf("  [%d] n_obs=%g features=%g  path=%s\n",
              i, flat[[i]]$n_obs, flat[[i]]$features,
              substr(flat[[i]]$path, 2, 80)))
}
cat("\n")

# ==============================================================================
# PART 5: CHECKPOINT -- filter out what is already done
# ==============================================================================

if (!dir.exists(cartella_robust)) dir.create(cartella_robust, recursive = TRUE)
checkpoint <- if (file.exists(checkpoint_file)) readLines(checkpoint_file) else character(0)

da_fare   <- Filter(function(x) !(x$path %in% checkpoint), flat)
gia_fatti <- length(flat) - length(da_fare)
cat("Already completed:", gia_fatti, "| To do:", length(da_fare), "\n\n")

# Split by size: large (n_obs * features >= 5000) vs small
da_fare_grandi  <- Filter(function(x) x$n_obs * x$features >= 5000, da_fare)
da_fare_piccoli <- Filter(function(x) x$n_obs * x$features <  5000, da_fare)
cat(sprintf("  Large (n*p >= 5000): %d | Small: %d\n\n",
            length(da_fare_grandi), length(da_fare_piccoli)))

# Write the total to the log (only if it doesn't exist yet)
if (!file.exists(log_file)) {
  writeLines(paste0("TOTALE: ", n_totali), log_file)
}

if (length(da_fare) == 0) {
  cat("All datasets already completed.\n")
  stop("Nothing to do -- go to PART 6 to aggregate.")
}


# ==============================================================================
# PART 6: SETUP (Linux: env vars only; Windows: makeCluster)
# ==============================================================================

data.table::setDTthreads(1)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

funzioni_sword <- c(
  "get_deviance", "select_top_n_SVM", "coef_unscaled_SVM",
  "best_split_SVM_ROT", "node_split_SVM_ROT", "tree_grow_SVM_ROT",
  "SVM.ROT", "Pred_coef_SVM", "precompute_predictions_SVM",
  "go_tree_svm_rot_optimized", "tree.predict_SpSVM_cor_new",
  "stripname", "depth_strip", "Depth_tree", "N_of_Leafs",
  "Leaf_or_Depth_for_CV", "depth_svm_rot", "toTree",
  "SVM.ROT.PRED.RF_all", "SVM.ROT.PRED.RF",
  "SVM.ROT.RF.OOB",
  "mean_variable_importance_SVM_ROT",
  "variable_importance_SVM_ROT",
  "Table_SVM.ROT", "Table_SVM.ROT_scaled", "coef_scaled",
  "plot_SVM.ROT"
)

funzioni_helper <- c(
  "correlazione_alberi_rf", "strength_rf_breiman",
  "RGA_metric", "parse_path_info", "%||%"
)

if (.Platform$OS.type == "windows") {
  library(doParallel)
  cl_large <- makeCluster(n_workers_large)
  .setup_cluster <- function(cl) {
    clusterExport(cl, varlist = "PROJ")
    clusterEvalQ(cl, {
      data.table::setDTthreads(1)
      Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
      library(data.table); library(infotheo); library(WeightSVM)
      library(data.tree);  library(stringr);  library(Metrics)
      library(fastDummies)
      `%||%` <- function(x, y) if (is.null(x)) y else x
      source(file.path(PROJ, "R", "SVM.ROT_Functions_ADAC.R"))
    })
    clusterExport(cl, varlist = c(funzioni_helper,
                                  "cartella_robust", "log_file", "checkpoint_file",
                                  "rf_var_fracs", "PARAM_RELATION",
                                  "PARAM_RAND_NTOPCOR", "PARAM_COST_C"))
  }
  .setup_cluster(cl_large)
  registerDoParallel(cl_large)
}


# ==============================================================================
# PART 7: PARALLELISM -- mclapply on Linux, foreach on Windows
#
# LINUX: mclapply with mc.preschedule = FALSE
#   Each dataset gets a FRESH fork of the parent process.
#   When the fork exits, the kernel frees ALL of its RAM.
#   No accumulation across successive datasets.
#
# Windows: foreach %dopar% with makeCluster (persistent workers)
# ==============================================================================

cat("========================================\n")
cat("STARTING SWORD SIMULATION (robust)\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

t_start <- Sys.time()

# --- Function that processes ONE dataset (used by both mclapply and foreach) ---
.process_item <- function(item) {
  gc()
  df   <- item$df
  path <- item$path

  # --- Output file path ---
  file_name <- paste0("robust__", gsub("\\$", "__", substr(path, 2, 999)))
  rds_path  <- file.path(cartella_robust, paste0(file_name, ".rds"))

  # Skip if it already exists (double protection on top of the checkpoint)
  if (file.exists(rds_path)) return(invisible(NULL))

  # --- Metadata ---
  info     <- parse_path_info(path)
  t0_total <- proc.time()

  # --- rf_var based on the ORIGINAL features (before dummy encoding) ---
  ncolumns   <- ncol(df) - 1
  var_per_rf <- unique(pmax(2, ceiling(ncolumns * rf_var_fracs)))

  # --- Dummy encoding if needed ---
  if (any(sapply(df, function(col) is.factor(col) || is.character(col)))) {
    df <- fastDummies::dummy_cols(
      df, remove_selected_columns = TRUE, remove_first_dummy = TRUE
    )
  }

  param_grid <- expand.grid(
    rand_ntopcor = PARAM_RAND_NTOPCOR,
    relation     = PARAM_RELATION,
    Cost_C       = PARAM_COST_C,
    rf_var       = var_per_rf,
    stringsAsFactors = FALSE
  )
  # Typically 36 rows (2x2x3x3); fewer if ncolumns is small

  # --- Serial loop over the grid ---
  grid_ris <- vector("list", nrow(param_grid))

  for (i in seq_len(nrow(param_grid))) {
    p  <- param_grid[i, ]
    ti <- proc.time()

    mod <- tryCatch(
      SVM.ROT.RF.OOB(
        Covariates    = df[, -1, drop = FALSE],
        y             = df[, 1],
        nmin          = 5,
        cp            = 0.00,
        n_perc        = 1,
        n_topCor      = 2,
        threshold_COR = 1,
        m             = 100,
        rf_var        = p$rf_var,
        rand_ntopcor  = p$rand_ntopcor,
        relation      = p$relation,
        Weight_Scheme = "robust",
        type_of_svm   = "C-classification",
        cost_C        = p$Cost_C,
        cost_nu       = 0.5,
        seed_BS       = 25,
        parallel      = FALSE,
        chunk         = FALSE,
        n_chunks      = NULL,
        OOB           = TRUE,
        n_workers     = NULL
      ),
      error = function(e) {
        cat(file = log_file, append = TRUE,
            sprintf("%s | ERROR | %s | combo %d | %s\n",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    file_name, i, e$message))
        NULL
      }
    )

    elapsed_i <- as.numeric((proc.time() - ti)[["elapsed"]])

    if (!is.null(mod)) {
      cor_tree <- correlazione_alberi_rf(mod$OOB_matrix)
      str_r2   <- strength_rf_breiman(mod$OOB_matrix, df[, 1])

      # RGA -- average OOB predictions
      y_oob   <- rowMeans(mod$OOB_matrix, na.rm = TRUE)
      rga_val <- RGA_metric(df[, 1], y_oob)

      res_row <- cbind(
        data.frame(
          relation     = p$relation,
          rand_ntopcor = p$rand_ntopcor,
          rf_var       = p$rf_var,
          W_scheme     = "robust",
          Cost_C       = p$Cost_C,
          sd_df        = sd(df[, 1], na.rm = TRUE),
          mse          = mod$MSE      %||% NA_real_,
          rmse         = mod$RMSE     %||% NA_real_,
          mae          = mod$MAE      %||% NA_real_,
          rsquared     = mod$Rsquared %||% NA_real_,
          rga          = rga_val,
          cor          = cor_tree,
          str_r2       = str_r2,
          time_grid_sec = elapsed_i,
          stringsAsFactors = FALSE
        ),
        info
      )

      # Tree and OOB timing
      if (!is.null(mod$time_tree)) {
        tt <- mod$time_tree
        res_row$time_tree_mean <- mean(tt, na.rm = TRUE)
        res_row$time_tree_sd   <- sd(tt,   na.rm = TRUE)
        res_row$time_tree_sum  <- sum(tt,   na.rm = TRUE)
      }
      if (!is.null(mod$time_oob)) {
        to <- mod$time_oob
        res_row$time_oob_mean <- mean(to, na.rm = TRUE)
        res_row$time_oob_sd   <- sd(to,   na.rm = TRUE)
        res_row$time_oob_sum  <- sum(to,   na.rm = TRUE)
      }

      # VI -- seed_1 only
      vi_norm <- if (!is.na(info$seed) && info$seed == 1) {
        vi_raw <- tryCatch(mean_variable_importance_SVM_ROT(mod), error = function(e) NULL)
        if (!is.null(vi_raw) && sum(vi_raw, na.rm = TRUE) > 0)
          vi_raw / sum(vi_raw, na.rm = TRUE)
        else NULL
      } else NULL

      grid_ris[[i]] <- list(res        = res_row,
                            OOB_matrix = mod$OOB_matrix,
                            time_tree  = mod$time_tree,
                            time_oob   = mod$time_oob,
                            VI         = vi_norm)
      # Free the forest right away (heavy): OOB_matrix is already copied into grid_ris
      rm(mod); gc()
    }
  } # end grid loop

  tempo_tot <- as.numeric((proc.time() - t0_total)[["elapsed"]])

  # --- Save .rds and checkpoint ---
  # y saved once outside the combos (identical across all 36)
  saveRDS(list(y             = df[, 1],
               combos        = grid_ris,
               tempo_tot_sec = tempo_tot),
          file = rds_path)
  write(path, file = checkpoint_file, append = TRUE)

  cat(file = log_file, append = TRUE,
      sprintf("%s | DONE | %s | %.1f min\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              file_name, tempo_tot / 60))

  invisible(NULL)
} # end .process_item


# --- Parallel execution ---
if (.Platform$OS.type != "windows") {

  # LINUX: mclapply with mc.preschedule = FALSE
  # Each dataset -> fresh fork -> RAM reclaimed when the fork exits
  if (length(da_fare_grandi) > 0) {
    cat(sprintf("\n--- Pass 1: %d large datasets (%d workers) ---\n",
                length(da_fare_grandi), n_workers_large))
    parallel::mclapply(da_fare_grandi, .process_item,
                       mc.cores       = n_workers_large,
                       mc.preschedule = FALSE)
  }
  if (length(da_fare_piccoli) > 0) {
    cat(sprintf("\n--- Pass 2: %d small datasets (%d workers) ---\n",
                length(da_fare_piccoli), n_workers_small))
    parallel::mclapply(da_fare_piccoli, .process_item,
                       mc.cores       = n_workers_small,
                       mc.preschedule = FALSE)
  }

} else {

  # WINDOWS: foreach %dopar% with makeCluster (persistent workers)
  .run_foreach_win <- function(lista_items, cl_win) {
    if (length(lista_items) == 0) return(invisible(NULL))
    registerDoParallel(cl_win)
    foreach(
      idx = seq_along(lista_items),
      .packages = c("fastDummies", "data.table", "infotheo",
                    "WeightSVM", "data.tree", "stringr", "Metrics")
    ) %dopar% { .process_item(lista_items[[idx]]) }
  }

  if (length(da_fare_grandi) > 0) {
    cat(sprintf("\n--- Pass 1: %d large datasets (%d workers) ---\n",
                length(da_fare_grandi), n_workers_large))
    .run_foreach_win(da_fare_grandi, cl_large)
  }
  if (length(da_fare_piccoli) > 0) {
    cat(sprintf("\n--- Pass 2: %d small datasets (%d workers) ---\n",
                length(da_fare_piccoli), n_workers_small))
    cl_small <- makeCluster(n_workers_small)
    .setup_cluster(cl_small)
    .run_foreach_win(da_fare_piccoli, cl_small)
    parallel::stopCluster(cl_small)
  }
  parallel::stopCluster(cl_large)

}

t_end <- Sys.time()
cat("\n========================================\n")
cat("SIMULATION COMPLETE!\n")
cat("Total time:", round(as.numeric(difftime(t_end, t_start, units = "hours")), 2), "h\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")


# ==============================================================================
# PART 8: AGGREGATION -- df_finale_SWORD_robust
# ==============================================================================

cat("Aggregating results into a data frame...\n")

file_rds <- list.files(cartella_robust, pattern = "\\.rds$", full.names = TRUE)
file_rds <- file_rds[basename(file_rds) != "risultati_completi.rds"]

df_finale_SWORD_robust <- do.call(rbind, lapply(file_rds, function(f) {
  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(obj)) return(NULL)
  # New structure: list(y, combos, tempo_tot_sec)
  combo_list <- if (!is.null(obj$combos)) obj$combos else obj[sapply(obj, is.list)]
  if (length(combo_list) == 0) return(NULL)
  do.call(rbind, lapply(combo_list, function(x) if (!is.null(x$res)) x$res else NULL))
}))

cat("Rows in df_finale_SWORD_robust:", nrow(df_finale_SWORD_robust), "\n")

save(df_finale_SWORD_robust,
     file = file.path(cartella_robust, "df_finale_SWORD_robust.rda"))
cat("Saved: risultati_full_simulation_robust/df_finale_SWORD_robust.rda\n")
cat("\nDONE.\n")
