# ==============================================================================
# FULL SIMULATION — Random Forest baseline, all seeds
# Parallelized over DATASETS, serial trees inside each worker
# Real-time progress log: use check_progress() from another console
#
# NOTE: this is the "from scratch" generation script (slow, meant for a
# many-core machine), included for full transparency/provenance. The light
# caches actually shipped (risultati_RF_simulation/df_finale_RF_ADAC.rda and
# VI_RF_seed1_cached.rds) are what the paper's figures/tables read from --
# this script is NOT required to reproduce any of them. Run this only to
# regenerate the RF baseline from scratch.
# ==============================================================================
rm(list = ls())
gc()
# ==============================================================================
# PART 0: SETUP
# ==============================================================================

# --- Install missing packages (first run only) ---
pkgs_needed <- c("data.table", "infotheo", "WeightSVM", "data.tree",
                 "stringr", "Metrics", "fastDummies",
                 "future", "furrr", "progressr",
                 "foreach", "doParallel",
                 "dplyr", "tidyr", "ggplot2", "ggpubr", "latex2exp",
                 "randomForest")

pkgs_missing <- pkgs_needed[!sapply(pkgs_needed, requireNamespace, quietly = TRUE)]
if (length(pkgs_missing) > 0) {
  cat("Installing missing packages:", paste(pkgs_missing, collapse = ", "), "\n")
  install.packages(pkgs_missing, repos = "https://cloud.r-project.org")
}

# --- Load packages ---
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

# --- Source functions ---
source(file.path(PROJ, "R", "SVM.ROT_Functions_ADAC.R"))


# --- Load simulated datasets (all seeds) ---
cat("Loading simulated datasets...\n")
load(file.path(DATA, "simulated_datasets_ADAC.rda"))
#list_simulated_df <- list_simulated_df["seed_9"]


# ==============================================================================
# PART 1: check_progress() FUNCTION
# Call it from the RStudio console while the simulation is running
# ==============================================================================

check_progress <- function(cartella = file.path(SCRIPT_DIR, "risultati_RF_simulation")) {

  # Count completed .rds files
  rds <- list.files(cartella, pattern = "\\.rds$")
  rds <- rds[rds != "risultati_completi.rds"]
  n_done <- length(rds)

  # Read the log
  log_file <- file.path(cartella, "progress_log.txt")
  if (file.exists(log_file)) {
    log_lines <- readLines(log_file, warn = FALSE)
    n_log <- length(log_lines)

    # Read n_totali from the first line
    n_totali <- NA
    if (n_log > 0 && grepl("^TOTALE:", log_lines[1])) {
      n_totali <- as.integer(sub("TOTALE: ", "", log_lines[1]))
      log_lines <- log_lines[-1]  # drop header
    }

    cat("============================================\n")
    if (!is.na(n_totali)) {
      pct <- round(100 * n_done / n_totali, 1)
      cat(sprintf("  Progress: %d / %d datasets (%.1f%%)\n", n_done, n_totali, pct))

      # Estimate remaining time
      done_lines <- log_lines[grepl("^\\d{4}-", log_lines)]
      if (length(done_lines) >= 2) {
        # Take timestamp of first and last completed
        first_ts <- as.POSIXct(substr(done_lines[1], 1, 19), format = "%Y-%m-%d %H:%M:%S")
        last_ts  <- as.POSIXct(substr(tail(done_lines, 1), 1, 19), format = "%Y-%m-%d %H:%M:%S")
        elapsed_min <- as.numeric(difftime(last_ts, first_ts, units = "mins"))
        if (elapsed_min > 0 && n_done > 0) {
          rate <- n_done / elapsed_min  # datasets/min
          remaining <- (n_totali - n_done) / rate
          cat(sprintf("  Rate: %.1f datasets/min\n", rate))
          cat(sprintf("  Estimated remaining time: %.0f min (%.1f h)\n", remaining, remaining / 60))
        }
      }
    } else {
      cat(sprintf("  Datasets completed: %d\n", n_done))
    }
    cat("============================================\n")

    # Last 10 completed
    done_lines <- log_lines[grepl("^\\d{4}-", log_lines)]
    if (length(done_lines) > 0) {
      cat("\nLast completed:\n")
      cat(paste(" ", tail(done_lines, 10), collapse = "\n"), "\n")
    }
  } else {
    cat("Simulation not started yet or empty folder.\n")
    cat("Found .rds datasets:", n_done, "\n")
  }

  invisible(n_done)
}

# helper
`%||%` <- function(x, y) if (is.null(x)) y else x

correlazione_alberi_rf <- function(OOB_matrix) {
  m <- ncol(OOB_matrix)
  cor_matrix <- matrix(NA, m, m)
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      valid_idx <- !is.na(OOB_matrix[, i]) & !is.na(OOB_matrix[, j])
      if (sum(valid_idx) > 2) {
        x <- OOB_matrix[valid_idx, i]
        y <- OOB_matrix[valid_idx, j]
        if (sd(x) > 0 && sd(y) > 0) {
          cor_val <- cor(x, y, use = "pairwise.complete.obs")
          cor_matrix[i, j] <- cor_val
          cor_matrix[j, i] <- cor_val
        }
      }
    }
  }
  diag(cor_matrix) <- 1
  mean(cor_matrix[upper.tri(cor_matrix)], na.rm = TRUE)
}

strength_rf_breiman <- function(OOB_matrix, target) {
  T <- ncol(OOB_matrix)
  pred_rf <- rowMeans(OOB_matrix, na.rm = TRUE)
  mse_rf <- cor(target, pred_rf, use = "pairwise.complete.obs")^2
  mse_trees <- sapply(1:T, function(t) {
    valid_idx <- !is.na(OOB_matrix[, t])
    if (any(valid_idx)) cor(target[valid_idx], OOB_matrix[valid_idx, t], use = "pairwise.complete.obs")^2 else NA
  })
  mean(mse_trees, na.rm = TRUE)
}


grid_train_RF_base <- function(df) {
  library(randomForest)
  library(Metrics)

  if (any(sapply(df, function(col) is.factor(col) || is.character(col)))) {
    df <- fastDummies::dummy_cols(df, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
  }

  set.seed(25)
  start_time <- Sys.time()
  rf <- randomForest(target ~ ., data = df, ntree = 100, keep.inbag = TRUE)
  end_time <- Sys.time()

  elapsed_min <- as.numeric(difftime(end_time, start_time, units = "mins"))
  y <- df$target
  pred_rf <- rf$predicted

  oob_predictions <- predict(rf, df, predict.all = TRUE, type = "response")
  preds_all <- oob_predictions$individual
  inbag <- rf$inbag
  n <- nrow(df)
  T <- ncol(preds_all)
  OOB_matrix <- matrix(NA, nrow = n, ncol = T)
  for (t in 1:T) {
    oob_idx <- inbag[, t] == 0
    OOB_matrix[oob_idx, t] <- preds_all[oob_idx, t]
  }

  mse <- mean((y - pred_rf)^2)
  rmse <- sqrt(mse)
  rsq <- (cor(y, pred_rf))^2
  mad_val <- mae(y, pred_rf)

  cor_tree <- correlazione_alberi_rf(OOB_matrix)
  strength <- strength_rf_breiman(OOB_matrix, df$target)

  list(
    res = data.frame(
      mse = mse, rmse = rmse, mae = mad_val, rsquared = rsq,
      cor = cor_tree, str_r2 = strength, sd_df = sd(y),
      time_min = elapsed_min, stringsAsFactors = FALSE
    ),
    mse_by_tree = rf$mse,
    VI = rf$importance / sum(rf$importance),
    OOB_matrix = OOB_matrix
  )
}


# ==============================================================================
# applica_modello_RF_GRID
# ==============================================================================
applica_modello_RF_GRID <- function(lista, path = "", cartella = file.path(SCRIPT_DIR, "risultati_RF_simulation")) {
  if (!dir.exists(cartella)) dir.create(cartella)
  checkpoint_file <- file.path(cartella, "checkpoint.txt")
  checkpoint <- if (file.exists(checkpoint_file)) readLines(checkpoint_file) else character(0)

  if (is.data.frame(lista)) {
    file_name <- gsub("\\$", "__", substr(path, 2, 999))
    rds_path <- file.path(cartella, paste0(file_name, ".rds"))

    if (path %in% checkpoint && file.exists(rds_path)) {
      cat("Already done:", substr(path, 2, 999), "\n")
      return(readRDS(rds_path))
    }

    elementi <- strsplit(path, "\\$")[[1]]
    elementi <- elementi[nzchar(elementi)]
    valori <- rep(NA, 7)
    if (length(elementi) == 7) {
      valori <- as.numeric(gsub("[^0-9.]", "", elementi))
    }
    info <- data.frame(
      seed = valori[1], n_obs = valori[2], features = valori[3],
      noise = valori[4], nonlin = valori[5], cat = valori[6],
      error_sd = valori[7]
    )

    ris <- tryCatch({
      risultato <- grid_train_RF_base(lista)
      risultato$res <- cbind(risultato$res, info)
      saveRDS(list(risultato), rds_path)
      write(path, file = checkpoint_file, append = TRUE)
      cat("Saved:", substr(path, 2, 999), "\n")
      list(risultato)
    }, error = function(e) {
      warning("Error in ", substr(path, 2, 999), ": ", e$message)
      return(NULL)
    })

    return(ris)
  } else if (is.list(lista)) {
    result_list <- list()
    for (n in names(lista)) {
      new_path <- paste0(path, "$", n)
      if (is.atomic(lista[[n]])) next
      result_list[[n]] <- applica_modello_RF_GRID(lista[[n]], new_path, cartella)
    }
    return(result_list)
  } else {
    return(NULL)
  }
}


# ==============================================================================
# PART 2: RF SIMULATION
# ==============================================================================

cat("========================================\n")
cat("STARTING RF SIMULATION\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

library(randomForest)

cartella_rf <- file.path(SCRIPT_DIR, "risultati_RF_simulation")

t_start_rf <- Sys.time()

risultati_rf <- applica_modello_RF_GRID(
  list_simulated_df,
  cartella = cartella_rf
)

t_end_rf <- Sys.time()
cat("RF completed in:", round(as.numeric(difftime(t_end_rf, t_start_rf, units = "mins")), 1), "min\n\n")

# Aggregate RF into a data frame
file_rds_rf <- list.files(cartella_rf, pattern = "\\.rds$", full.names = TRUE)
df_finale_RF <- do.call(rbind, lapply(file_rds_rf, function(f) {
  res <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(res)) return(NULL)
  do.call(rbind, lapply(res, function(x) x$res))
}))

cat("Rows in df_finale_RF:", nrow(df_finale_RF), "\n")

save(df_finale_RF, file = file.path(cartella_rf, "df_finale_RF_ADAC.rda"))
cat("Saved: df_finale_RF_ADAC.rda\n\n")
