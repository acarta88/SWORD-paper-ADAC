# ==============================================================================
# Table 2 (tab:stability_bootstrap) — OIW-VI STABILITY ACROSS TREES: SWORD vs RF
# Per-tree variable-importance stability under correlated features
# ------------------------------------------------------------------------------
# Design: simulated dataset with n=200, p=10 (5 informative, 5 noise)
# Two within-group correlation levels: rho = 0 (baseline) and rho = 0.95
# For each tree b=1:100:
#   - SWORD: OIW-VI_bj computed on all internal nodes via
#            variable_importance_SVM_ROT()  (pre_splitDev - dev)
#   - RF:    MDI_bj computed on a single bootstrap tree via randomForest
# Stability measure: CV = SD/mean of the normalized VI across the B trees
# ==============================================================================

library(MASS)
library(randomForest)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

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
dir_out <- file.path(PROJ, "tables")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)
if (!interactive()) pdf(NULL)

source(file.path(PROJ, "R", "SVM.ROT_Functions_ADAC.R"))

set.seed(42)


# ==============================================================================
# 1. PARAMETERS
# ==============================================================================

n        <- 200
p_info   <- 5
p_noise  <- 5
p        <- p_info + p_noise
B        <- 100
rho_vals <- c(0, 0.95)


# ==============================================================================
# 2. FUNCTION: per-tree OIW-VI SWORD (all nodes, normalized)
# ==============================================================================

vi_sword_per_tree <- function(foresta) {
  do.call(rbind, Filter(Negate(is.null), lapply(seq_along(foresta$trees), function(b) {
    vi <- tryCatch(
      variable_importance_SVM_ROT(foresta$trees[[b]]),
      error = function(e) NULL
    )
    if (is.null(vi) || sum(vi, na.rm = TRUE) == 0) return(NULL)
    vi_norm <- vi / sum(vi, na.rm = TRUE)
    data.frame(albero = b, variable = names(vi_norm), VI = as.numeric(vi_norm),
               row.names = NULL)
  })))
}


# ==============================================================================
# 3. FUNCTION: per-tree MDI RF (single bootstrap trees, normalized)
# ==============================================================================

vi_rf_per_tree <- function(df, B) {
  n <- nrow(df)
  do.call(rbind, Filter(Negate(is.null), lapply(seq_len(B), function(b) {
    boot_idx <- sample(n, n, replace = TRUE)
    rf_b <- tryCatch(
      randomForest(target ~ ., data = df[boot_idx, ], ntree = 1,
                   importance = FALSE, keep.forest = FALSE),
      error = function(e) NULL
    )
    if (is.null(rf_b)) return(NULL)
    vi <- importance(rf_b, type = 2)[, 1]
    tot <- sum(vi, na.rm = TRUE)
    if (tot == 0) return(NULL)
    data.frame(albero = b, variable = names(vi), VI = as.numeric(vi / tot),
               row.names = NULL)
  })))
}


# ==============================================================================
# 4. LOOP OVER CORRELATION LEVELS
# ==============================================================================
set.seed(10)
risultati <- list()

for (rho in rho_vals) {
  cat(sprintf("\n=== rho = %.2f ===\n", rho))

  # Simulate data
  Sigma_info <- matrix(rho, p_info, p_info); diag(Sigma_info) <- 1
  X_info  <- mvrnorm(n, mu = rep(0, p_info), Sigma = Sigma_info)
  X_noise <- matrix(rnorm(n * p_noise), nrow = n)
  X       <- cbind(X_info, X_noise)
  colnames(X) <- c(paste0("x_info", 1:p_info), paste0("x_noise", 1:p_noise))
  y  <- X_info %*% rep(1, p_info) + rnorm(n)
  df <- data.frame(target = as.numeric(y), X)

  # SWORD
  cat("  Fitting SWORD...\n")
  foresta <- tryCatch(
    SVM.ROT.RF.OOB(
      Covariates    = df[, -1],
      y             = df[, 1],
      nmin          = 5,
      cp            = 0.00,
      m             = B,
      rf_var        = p,
      rand_ntopcor  = TRUE,
      relation      = "Pearson",
      Weight_Scheme = "scale",
      cost_C        = 1,
      seed_BS       = 25,
      OOB           = FALSE,
      parallel      = FALSE
    ),
    error = function(e) { cat("  ERROR SWORD:", e$message, "\n"); NULL }
  )

  if (!is.null(foresta)) {
    vi_sw <- vi_sword_per_tree(foresta)
    vi_sw$metodo <- "SWORD"
    vi_sw$rho    <- rho
    cat(sprintf("  SWORD: %d per-tree VI rows\n", nrow(vi_sw)))
  } else {
    vi_sw <- NULL
  }

  # RF
  cat("  Fitting RF (", B, "bootstrap)...\n")
  vi_rf <- vi_rf_per_tree(df, B)
  vi_rf$metodo <- "RF"
  vi_rf$rho    <- rho
  cat(sprintf("  RF: %d per-tree VI rows\n", nrow(vi_rf)))

  risultati[[as.character(rho)]] <- rbind(vi_sw, vi_rf)
}

df_vi <- do.call(rbind, risultati)
df_vi$var_type <- ifelse(grepl("info", df_vi$variable), "Informative", "Noise")
df_vi$rho_lab  <- paste0("rho == ", df_vi$rho)


# ==============================================================================
# 5. CV PER VARIABLE
# ==============================================================================

cv_per_var <- df_vi %>%
  group_by(rho, metodo, variable, var_type) %>%
  summarise(
    vi_mean = mean(VI, na.rm = TRUE),
    vi_sd   = sd(VI,   na.rm = TRUE),
    cv      = ifelse(vi_mean > 1e-8, vi_sd / vi_mean, NA_real_),
    .groups = "drop"
  )

cat("\n=== CV per (rho, method, variable) ===\n")
print(cv_per_var %>% arrange(rho, metodo, var_type, variable), n = 80)

cat("\n=== Mean CV summary per (rho, method, var_type) ===\n")
print(
  cv_per_var %>%
    group_by(rho, metodo, var_type) %>%
    summarise(CV_mean   = round(mean(cv,   na.rm = TRUE), 3),
              CV_median = round(median(cv, na.rm = TRUE), 3),
              .groups = "drop") %>%
    arrange(rho, metodo, var_type)
)


# ==============================================================================
# 6. RANK STABILITY: Spearman(VI^(b), VI^(b')) between pairs of trees
# ==============================================================================

# For each (rho, method): B x p matrix of per-tree VI
# Compute Spearman correlation between all pairs of trees (b, b')

spearman_pairs <- do.call(rbind, lapply(rho_vals, function(r) {
  do.call(rbind, lapply(c("SWORD", "RF"), function(met) {

    # B x p matrix: each row = a tree's VI vector
    mat <- df_vi %>%
      filter(rho == r, metodo == met) %>%
      dplyr::select(albero, variable, VI) %>%
      pivot_wider(names_from = variable, values_from = VI) %>%
      arrange(albero) %>%
      dplyr::select(-albero) %>%
      as.matrix()

    alberi <- seq_len(nrow(mat))
    coppie <- combn(alberi, 2)   # 4950 pairs for B=100

    rho_vals_sp <- apply(coppie, 2, function(idx) {
      cor(mat[idx[1], ], mat[idx[2], ], method = "spearman", use = "complete.obs")
    })

    data.frame(
      rho    = r,
      metodo = met,
      sp_mean   = mean(rho_vals_sp,   na.rm = TRUE),
      sp_median = median(rho_vals_sp, na.rm = TRUE),
      sp_sd     = sd(rho_vals_sp,     na.rm = TRUE),
      sp_q10    = quantile(rho_vals_sp, 0.10, na.rm = TRUE),
      sp_values = I(list(rho_vals_sp))
    )
  }))
}))

cat("\n=== Rank stability: Spearman between pairs of trees ===\n")
print(spearman_pairs[, c("rho","metodo","sp_mean","sp_median","sp_sd","sp_q10")])

# Wilcoxon SWORD vs RF for each rho
for (r in rho_vals) {
  sw <- unlist(spearman_pairs$sp_values[spearman_pairs$rho == r & spearman_pairs$metodo == "SWORD"])
  rf <- unlist(spearman_pairs$sp_values[spearman_pairs$rho == r & spearman_pairs$metodo == "RF"])
  wt <- wilcox.test(sw, rf)
  cat(sprintf("rho=%.2f: Average Spearman SWORD=%.3f vs RF=%.3f  (Wilcoxon W=%.0f, p=%.2e)\n",
              r, mean(sw), mean(rf), wt$statistic, wt$p.value))
}

