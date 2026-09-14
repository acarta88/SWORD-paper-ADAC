# =============================================================================
# Table 1 (tab:wilcox) — Wilcoxon signed-rank tests on matched pairs — NRMSE
# (72 combinations); 243 x 72 matrix aggregated over seeds
#
# W_scheme: "robust" (risultati_full_simulation_robust)
#           "scale"      (risultati_full_simulation_scale)
#
# Effect of the 5 SWORD parameters on NRMSE (matched-pair):
#   rand_label (k=1), relation (k=1), Cost_C (k=3),
#   rf_var_frac (k=3), W_scheme (k=1)
#
# Output CSV:
#   tables/Wilcoxon_paired_results_72.csv
# =============================================================================

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
SCRIPT_DIR <- .get_script_dir()
PROJ    <- normalizePath(file.path(SCRIPT_DIR, ".."))
CACHE   <- file.path(PROJ, "results_cache")
dir_out <- file.path(PROJ, "tables")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)

# =============================================================================
# 1. Loading function
# =============================================================================

ricostruisci_df_new <- function(cartella) {
  file_rds <- list.files(cartella, pattern = "\\.rds$", full.names = TRUE)
  file_rds <- file_rds[basename(file_rds) != "risultati_completi.rds"]
  cat("  .rds files in", basename(cartella), ":", length(file_rds), "\n")
  do.call(rbind, Filter(Negate(is.null), lapply(file_rds, function(f) {
    obj <- tryCatch(readRDS(f), error = function(e) {
      cat("  ERROR:", basename(f), "\n"); NULL
    })
    if (is.null(obj)) return(NULL)
    combo_list <- if (!is.null(obj$combos)) obj$combos else
      Filter(function(x) is.list(x) && !is.null(x$res), obj)
    if (length(combo_list) == 0) return(NULL)
    do.call(rbind, Filter(Negate(is.null), lapply(combo_list, function(x) {
      if (is.null(x$res)) return(NULL)
      x$res
    })))
  })))
}

# =============================================================================
# 2. Data loading and preparation
# =============================================================================

path_cache <- file.path(CACHE, "df_SWORD72_cached.rds")

if (file.exists(path_cache)) {
  cat("Loading SWORD 72 (scale+robust) from cache...\n")
  df <- readRDS(path_cache)
} else {
  cat("Cache not found - rebuilding from .rds\n")

  cat("Loading SWORD robust...\n")
  df_rob <- ricostruisci_df_new("risultati_full_simulation_robust")
  if (!"W_scheme" %in% names(df_rob)) df_rob$W_scheme <- "robust"

  cat("Loading SWORD scale...\n")
  df_scl <- ricostruisci_df_new("risultati_full_simulation_scale")
  df_scl$W_scheme <- "scale"

  df <- bind_rows(df_rob, df_scl)
  cat("SWORD combined:", nrow(df), "rows\n")
  saveRDS(df, path_cache)
  cat("Cache saved to", path_cache, "\n")
}

df <- df %>%
  mutate(
    features    = as.integer(features),
    rf_var      = as.integer(rf_var),
    relation    = factor(relation, levels = c("Pearson", "MI")),
    Cost_C      = factor(as.character(Cost_C), levels = c("0.1", "1", "10")),
    rf_var_frac = factor(
      cut(rf_var / features, breaks = c(0, 0.45, 0.75, Inf),
          labels = c("30%", "60%", "100%"), right = TRUE),
      levels = c("30%", "60%", "100%")
    ),
    rand_label  = factor(ifelse(as.logical(rand_ntopcor), "Random", "Fixed"),
                         levels = c("Fixed", "Random")),
    W_scheme    = factor(W_scheme),
    NRMSE       = rmse / sd_df,
    dgp_id      = paste(n_obs, features, noise, nonlin, cat, error_sd, sep = "_")
  )

# Aggregation: mean over seeds -> 243 x 72 rows (long format)
df_agg <- df %>%
  group_by(dgp_id, W_scheme, relation, Cost_C, rf_var_frac, rand_label) %>%
  summarise(NRMSE_mean = mean(NRMSE, na.rm = TRUE),
            n_seed     = n(), .groups = "drop")

cat(sprintf("Aggregated rows: %d  (expected 243x72 = %d)\n",
            nrow(df_agg), 243 * 72))
cat(sprintf("Unique DGPs: %d | average seeds per cell: %.1f\n\n",
            n_distinct(df_agg$dgp_id), mean(df_agg$n_seed)))

# =============================================================================
# 3. Wilcoxon helper function
# =============================================================================

wilcox_paired <- function(x, y, label_x, label_y, bonf_k = 1) {
  diffs  <- x - y
  wt     <- wilcox.test(x, y, paired = TRUE, exact = FALSE, correct = TRUE)
  N      <- length(diffs)

  W_plus <- as.numeric(wt$statistic)
  mu_W   <- N * (N + 1) / 4
  sig_W  <- sqrt(N * (N + 1) * (2 * N + 1) / 24)
  Z      <- (W_plus - mu_W) / sig_W * sign(median(diffs))

  r_ros  <- min(abs(Z) / sqrt(N), 1)
  W_tot  <- N * (N + 1) / 2
  r_rb   <- (2 * W_plus - W_tot) / W_tot

  pct_x_wins <- 100 * mean(diffs > 0)
  p_raw  <- wt$p.value
  p_bon  <- min(p_raw * bonf_k, 1)

  data.frame(
    Confronto   = paste(label_x, "vs", label_y),
    Best        = ifelse(median(diffs) > 0, label_y, label_x),
    N_coppie    = N,
    pct_x_wins  = round(pct_x_wins, 1),
    median_diff = round(median(diffs), 5),
    r_Rosenthal = round(r_ros, 3),
    r_rb        = round(r_rb, 3),
    p_raw       = signif(p_raw, 3),
    p_bonf      = signif(p_bon, 3),
    sig         = ifelse(p_bon < 0.001, "***",
                  ifelse(p_bon < 0.01,  "**",
                  ifelse(p_bon < 0.05,  "*", "n.s.")))
  )
}

# =============================================================================
# 4. W_scheme: robust vs scale  (k=1)
# =============================================================================

cat("=== W_scheme: robust vs scale ===\n")
pw_ws <- df_agg %>%
  dplyr::select(dgp_id, relation, Cost_C, rf_var_frac, rand_label,
                W_scheme, NRMSE_mean) %>%
  pivot_wider(names_from = W_scheme, values_from = NRMSE_mean) %>%
  drop_na(robust, scale)

cat(sprintf("Valid pairs: %d\n", nrow(pw_ws)))
res_ws <- wilcox_paired(pw_ws$robust, pw_ws$scale,
                        "robust", "scale", bonf_k = 1)
print(res_ws)

# =============================================================================
# 5. rand_label: Fixed vs Random  (k=1)
# =============================================================================

cat("\n=== rand_label: Fixed vs Random ===\n")
pw_rand <- df_agg %>%
  dplyr::select(dgp_id, W_scheme, relation, Cost_C, rf_var_frac,
                rand_label, NRMSE_mean) %>%
  pivot_wider(names_from = rand_label, values_from = NRMSE_mean) %>%
  drop_na(Fixed, Random)

cat(sprintf("Valid pairs: %d\n", nrow(pw_rand)))
res_rand <- wilcox_paired(pw_rand$Fixed, pw_rand$Random,
                          "Fixed", "Random", bonf_k = 1)
print(res_rand)

# =============================================================================
# 6. relation: Pearson vs MI  (k=1)
# =============================================================================

cat("\n=== relation: Pearson vs MI ===\n")
pw_rel <- df_agg %>%
  dplyr::select(dgp_id, W_scheme, Cost_C, rf_var_frac, rand_label,
                relation, NRMSE_mean) %>%
  pivot_wider(names_from = relation, values_from = NRMSE_mean) %>%
  drop_na(Pearson, MI)

cat(sprintf("Valid pairs: %d\n", nrow(pw_rel)))
res_rel <- wilcox_paired(pw_rel$Pearson, pw_rel$MI,
                         "Pearson", "MI", bonf_k = 1)
print(res_rel)

# =============================================================================
# 7. Cost_C: pairwise comparisons + Bonferroni (k=3)
# =============================================================================

cat("\n=== Cost_C: pairwise comparisons (Bonferroni k=3) ===\n")
livelli_c <- c("0.1", "1", "10")

res_c <- bind_rows(lapply(combn(livelli_c, 2, simplify = FALSE), function(pair) {
  pw <- df_agg %>%
    filter(Cost_C %in% pair) %>%
    dplyr::select(dgp_id, W_scheme, relation, rf_var_frac, rand_label,
                  Cost_C, NRMSE_mean) %>%
    pivot_wider(names_from = Cost_C, values_from = NRMSE_mean) %>%
    drop_na(all_of(pair))
  wilcox_paired(pw[[pair[1]]], pw[[pair[2]]],
                paste0("C=", pair[1]), paste0("C=", pair[2]),
                bonf_k = 3)
}))
print(res_c)

# =============================================================================
# 8. rf_var_frac: pairwise comparisons + Bonferroni (k=3)
# =============================================================================

cat("\n=== rf_var_frac: pairwise comparisons (Bonferroni k=3) ===\n")
livelli_rf <- c("30%", "60%", "100%")

res_rf <- bind_rows(lapply(combn(livelli_rf, 2, simplify = FALSE), function(pair) {
  pw <- df_agg %>%
    filter(rf_var_frac %in% pair) %>%
    dplyr::select(dgp_id, W_scheme, relation, Cost_C, rand_label,
                  rf_var_frac, NRMSE_mean) %>%
    pivot_wider(names_from = rf_var_frac, values_from = NRMSE_mean) %>%
    drop_na(all_of(pair))
  wilcox_paired(pw[[pair[1]]], pw[[pair[2]]],
                pair[1], pair[2], bonf_k = 3)
}))
print(res_rf)

# =============================================================================
# 9. Summary table (5 parameters)
# =============================================================================

cat("\n=== SUMMARY TABLE (5 parameters) ===\n")
res_all <- bind_rows(
  res_ws   %>% mutate(Parametro = "W_scheme"),
  res_rand %>% mutate(Parametro = "rand_label"),
  res_rel  %>% mutate(Parametro = "relation"),
  res_c    %>% mutate(Parametro = "Cost_C"),
  res_rf   %>% mutate(Parametro = "rf_var_frac")
) %>%
  dplyr::select(Parametro, Confronto, Best, N_coppie, pct_x_wins,
                median_diff, r_Rosenthal, r_rb, p_raw, p_bonf, sig)

print(res_all, row.names = FALSE)

# =============================================================================
# 11. Save CSV
# =============================================================================

write.csv(res_all,
          file.path(dir_out, "Wilcoxon_paired_results_72.csv"),
          row.names = FALSE)


