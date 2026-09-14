# =============================================================================
# Table 3 (tab:wilcox_corstr_params) + Table A1 (tab:wilcox_72vsRF)
# Wilcoxon signed-rank tests — cor and str_r2  (72 combinations)
#
# Part A -> Table 3: Effect of the 5 SWORD parameters on cor and str_r2 (matched-pair)
#          Parameters: rand_label, relation, Cost_C, rf_var_frac, W_scheme
#
# Part B -> Table A1: 72 SWORD vs RF configurations (matched-pair on 243 DGP)
#          Bonferroni k=72, effect size: r Rosenthal + rank-biserial r_rb
#
# W_scheme: "robust" (original Weight_Scheme) | "scale" (Weight_Scheme=scale)
#
# Output CSV:
#   tables/Wilcoxon_corstr_5parametri_72.csv
#   tables/Wilcoxon_corstr_72vsRF.csv
#   tables/tab_wilcox_corstr_params_72.csv   (paper table tab:wilcox_corstr_params)
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
RF_DIR  <- file.path(SCRIPT_DIR, "risultati_RF_simulation")
dir_out <- file.path(PROJ, "tables")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)

# =============================================================================
# 1. Loading functions (two different structures)
# =============================================================================


ricostruisci_df_old <- function(cartella) {
  file_rds <- list.files(cartella, pattern = "\\.rds$", full.names = TRUE)
  file_rds <- file_rds[basename(file_rds) != "risultati_completi.rds"]
  cat("  .rds files in", basename(cartella), ":", length(file_rds), "\n")
  do.call(rbind, Filter(Negate(is.null), lapply(file_rds, function(f) {
    res <- tryCatch(readRDS(f), error = function(e) {
      cat("  ERROR:", basename(f), "\n"); NULL
    })
    if (is.null(res)) return(NULL)
    do.call(rbind, Filter(Negate(is.null), lapply(res, function(x) {
      if (is.null(x) || is.null(x$res)) return(NULL)
      x$res
    })))
  })))
}

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
# 2. Data loading
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








rf_cache_path <- file.path(RF_DIR, "df_finale_RF_ADAC.rda")
if (file.exists(rf_cache_path)) {
  cat("Loading RF from cache...\n")
  load(rf_cache_path)  # -> df_finale_RF
  df_rf <- df_finale_RF
} else {
  cat("RF cache not found - rebuilding from .rds (needs the per-dataset raw files)\n")
  df_rf <- ricostruisci_df_old(RF_DIR)
}

# --- Seeds common to all three ---
seeds_comuni <- Reduce(intersect, list(
  unique(df$seed),
  unique(df_rf$seed)
))
cat(sprintf("Common seeds: %s  (n=%d)\n",
            paste(sort(seeds_comuni), collapse = ", "), length(seeds_comuni)))

df    <- df[df$seed    %in% seeds_comuni, ]
df_rf <- df_rf[df_rf$seed %in% seeds_comuni, ]

# --- Derived SWORD variables ---
df <- df %>%
  mutate(
    features    = as.integer(features),
    rf_var      = as.integer(rf_var),
    relation    = factor(relation,    levels = c("Pearson", "MI")),
    Cost_C      = factor(as.character(Cost_C), levels = c("0.1", "1", "10")),
    rf_var_frac = factor(
      cut(rf_var / features, breaks = c(0, 0.45, 0.75, Inf),
          labels = c("30%", "60%", "100%"), right = TRUE),
      levels = c("30%", "60%", "100%")
    ),
    rand_label  = factor(ifelse(as.logical(rand_ntopcor), "Random", "Fixed"),
                         levels = c("Fixed", "Random")),
    W_scheme    = factor(W_scheme, levels = c("robust", "scale")),
    dgp_id      = paste(n_obs, features, noise, nonlin, cat, error_sd, sep = "_")
  )

# --- dgp_id for RF ---
df_rf <- df_rf %>%
  mutate(
    features = as.integer(features),
    dgp_id   = paste(n_obs, features, noise, nonlin, cat, error_sd, sep = "_")
  )

cat(sprintf("SWORD: %d rows | RF: %d rows\n", nrow(df), nrow(df_rf)))
cat(sprintf("SWORD — seed: %s | W_scheme: %s\n",
            paste(sort(unique(df$seed)),      collapse = ", "),
            paste(levels(df$W_scheme),        collapse = ", ")))
cat(sprintf("RF    — seed: %s\n\n",
            paste(sort(unique(df_rf$seed)), collapse = ", ")))

# =============================================================================
# 3. Aggregation over seeds: mean per DGP x configuration
# =============================================================================

# SWORD: 243 x 72 rows expected (243 DGP x 72 combos)
df_sw_agg <- df %>%
  group_by(dgp_id, W_scheme, relation, Cost_C, rf_var_frac, rand_label) %>%
  summarise(
    cor_mean    = mean(cor,    na.rm = TRUE),
    str_r2_mean = mean(str_r2, na.rm = TRUE),
    n_seed      = n(),
    .groups     = "drop"
  )

cat(sprintf("SWORD aggregated: %d rows (expected 243x72 = %d)\n",
            nrow(df_sw_agg), 243 * 72))

# RF: 243 rows expected
df_rf_agg <- df_rf %>%
  group_by(dgp_id) %>%
  summarise(
    cor_mean_rf    = mean(cor,    na.rm = TRUE),
    str_r2_mean_rf = mean(str_r2, na.rm = TRUE),
    n_seed         = n(),
    .groups        = "drop"
  )

cat(sprintf("RF aggregated:    %d rows (expected 243)\n\n", nrow(df_rf_agg)))

# =============================================================================
# 4. Wilcoxon helper function
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

  pct_x_higher <- 100 * mean(diffs > 0)
  p_raw  <- wt$p.value
  p_bon  <- min(p_raw * bonf_k, 1)

  data.frame(
    Confronto    = paste(label_x, "vs", label_y),
    Higher       = ifelse(median(diffs) > 0, label_x, label_y),
    N_coppie     = N,
    pct_x_higher = round(pct_x_higher, 1),
    median_diff  = round(median(diffs), 5),
    r_Rosenthal  = round(r_ros, 3),
    r_rb         = round(r_rb, 3),
    p_raw        = signif(p_raw, 3),
    p_bonf       = signif(p_bon, 3),
    sig          = ifelse(p_bon < 0.001, "***",
                   ifelse(p_bon < 0.01,  "**",
                   ifelse(p_bon < 0.05,  "*", "n.s.")))
  )
}

# =============================================================================
# PART A - EFFECT OF THE 5 PARAMETERS ON cor AND str_r2
#   parameters: rand_label (k=1), relation (k=1), Cost_C (k=3),
#              rf_var_frac (k=3), W_scheme (k=1)
# =============================================================================

run_5factor_tests <- function(df_agg, metric_col, metric_name) {

  sep <- paste(rep("=", 65), collapse = "")
  cat(sprintf("\n%s\nPART A — Effect of 5 parameters on %s\n%s\n\n",
              sep, metric_name, sep))

  # --- W_scheme: robust vs scale ---
  cat(sprintf("--- W_scheme: robust vs scale (%s, k=1) ---\n", metric_name))
  pw_ws <- df_agg %>%
    dplyr::select(dgp_id, relation, Cost_C, rf_var_frac, rand_label, W_scheme,
                  value = !!sym(metric_col)) %>%
    pivot_wider(names_from = W_scheme, values_from = value) %>%
    drop_na(robust, scale)
  cat(sprintf("Valid pairs: %d\n", nrow(pw_ws)))
  res_ws <- wilcox_paired(pw_ws$robust, pw_ws$scale,
                          "robust", "scale", bonf_k = 1) %>%
    mutate(Parametro = "W_scheme", Metrica = metric_name)
  print(res_ws)

  # --- rand_label: Fixed vs Random ---
  cat(sprintf("\n--- rand_label: Fixed vs Random (%s, k=1) ---\n", metric_name))
  pw_rand <- df_agg %>%
    dplyr::select(dgp_id, W_scheme, relation, Cost_C, rf_var_frac, rand_label,
                  value = !!sym(metric_col)) %>%
    pivot_wider(names_from = rand_label, values_from = value) %>%
    drop_na(Fixed, Random)
  cat(sprintf("Valid pairs: %d\n", nrow(pw_rand)))
  res_rand <- wilcox_paired(pw_rand$Fixed, pw_rand$Random,
                            "Fixed", "Random", bonf_k = 1) %>%
    mutate(Parametro = "rand_label", Metrica = metric_name)
  print(res_rand)

  # --- relation: Pearson vs MI ---
  cat(sprintf("\n--- relation: Pearson vs MI (%s, k=1) ---\n", metric_name))
  pw_rel <- df_agg %>%
    dplyr::select(dgp_id, W_scheme, Cost_C, rf_var_frac, rand_label, relation,
                  value = !!sym(metric_col)) %>%
    pivot_wider(names_from = relation, values_from = value) %>%
    drop_na(Pearson, MI)
  cat(sprintf("Valid pairs: %d\n", nrow(pw_rel)))
  res_rel <- wilcox_paired(pw_rel$Pearson, pw_rel$MI,
                           "Pearson", "MI", bonf_k = 1) %>%
    mutate(Parametro = "relation", Metrica = metric_name)
  print(res_rel)

  # --- Cost_C: pairwise Bonferroni k=3 ---
  cat(sprintf("\n--- Cost_C: pairwise Bonferroni k=3 (%s) ---\n", metric_name))
  livelli_c <- c("0.1", "1", "10")
  res_c <- bind_rows(lapply(combn(livelli_c, 2, simplify = FALSE), function(pair) {
    pw <- df_agg %>%
      filter(Cost_C %in% pair) %>%
      dplyr::select(dgp_id, W_scheme, relation, rf_var_frac, rand_label, Cost_C,
                    value = !!sym(metric_col)) %>%
      pivot_wider(names_from = Cost_C, values_from = value) %>%
      drop_na(all_of(pair))
    wilcox_paired(pw[[pair[1]]], pw[[pair[2]]],
                  paste0("C=", pair[1]), paste0("C=", pair[2]),
                  bonf_k = 3) %>%
      mutate(Parametro = "Cost_C", Metrica = metric_name)
  }))
  print(res_c)

  # --- rf_var_frac: pairwise Bonferroni k=3 ---
  cat(sprintf("\n--- rf_var_frac: pairwise Bonferroni k=3 (%s) ---\n", metric_name))
  livelli_rf <- c("30%", "60%", "100%")
  res_rf <- bind_rows(lapply(combn(livelli_rf, 2, simplify = FALSE), function(pair) {
    pw <- df_agg %>%
      filter(rf_var_frac %in% pair) %>%
      dplyr::select(dgp_id, W_scheme, relation, Cost_C, rand_label, rf_var_frac,
                    value = !!sym(metric_col)) %>%
      pivot_wider(names_from = rf_var_frac, values_from = value) %>%
      drop_na(all_of(pair))
    wilcox_paired(pw[[pair[1]]], pw[[pair[2]]],
                  pair[1], pair[2], bonf_k = 3) %>%
      mutate(Parametro = "rf_var_frac", Metrica = metric_name)
  }))
  print(res_rf)

  bind_rows(res_ws, res_rand, res_rel, res_c, res_rf)
}

res_A_cor <- run_5factor_tests(df_sw_agg, "cor_mean",    "cor")
res_A_str <- run_5factor_tests(df_sw_agg, "str_r2_mean", "str_r2")

res_A_all <- bind_rows(res_A_cor, res_A_str) %>%
  dplyr::select(Metrica, Parametro, Confronto, Higher, N_coppie, pct_x_higher,
                median_diff, r_Rosenthal, r_rb, p_raw, p_bonf, sig)

cat("\n\n=== SUMMARY TABLE PART A ===\n")
print(res_A_all, row.names = FALSE)

# =============================================================================
# PART B — 72 SWORD vs RF CONFIGURATIONS
# For each config: Wilcoxon paired on 243 DGP (Bonferroni k=72)
# =============================================================================

run_72vsRF <- function(df_sw_agg, df_rf_agg, sw_col, rf_col, metric_name) {

  sep <- paste(rep("=", 65), collapse = "")
  cat(sprintf("\n%s\nPART B — 72 SWORD vs RF configurations: %s  (Bonferroni k=72)\n%s\n\n",
              sep, metric_name, sep))

  configs <- df_sw_agg %>%
    distinct(W_scheme, relation, rand_label, Cost_C, rf_var_frac) %>%
    arrange(W_scheme, relation, rand_label, Cost_C, rf_var_frac)

  cat(sprintf("Configurations: %d\n\n", nrow(configs)))

  results <- lapply(seq_len(nrow(configs)), function(i) {
    cfg <- configs[i, ]

    sw_vals <- df_sw_agg %>%
      filter(W_scheme    == cfg$W_scheme,
             relation    == cfg$relation,
             rand_label  == cfg$rand_label,
             Cost_C      == cfg$Cost_C,
             rf_var_frac == cfg$rf_var_frac) %>%
      dplyr::select(dgp_id, value_sword = !!sym(sw_col))

    both <- sw_vals %>%
      inner_join(
        df_rf_agg %>% dplyr::select(dgp_id, value_rf = !!sym(rf_col)),
        by = "dgp_id"
      ) %>%
      drop_na(value_sword, value_rf)

    if (nrow(both) < 5) return(NULL)

    wilcox_paired(both$value_sword, both$value_rf,
                  "SWORD", "RF", bonf_k = 72) %>%
      mutate(
        Metrica     = metric_name,
        W_scheme    = as.character(cfg$W_scheme),
        relation    = as.character(cfg$relation),
        rand_label  = as.character(cfg$rand_label),
        Cost_C      = as.character(cfg$Cost_C),
        rf_var_frac = as.character(cfg$rf_var_frac),
        combo       = paste(cfg$W_scheme, cfg$relation, cfg$rand_label,
                            paste0("C=", cfg$Cost_C),
                            cfg$rf_var_frac, sep = " | ")
      )
  })

  bind_rows(results) %>%
    dplyr::select(Metrica, combo, W_scheme, relation, rand_label, Cost_C, rf_var_frac,
                  Higher, N_coppie, pct_x_higher, median_diff,
                  r_Rosenthal, r_rb, p_raw, p_bonf, sig)
}

res_B_cor <- run_72vsRF(df_sw_agg, df_rf_agg, "cor_mean",    "cor_mean_rf",    "cor")

res_B_str <- run_72vsRF(df_sw_agg, df_rf_agg, "str_r2_mean", "str_r2_mean_rf", "str_r2")

# Print results
cat("\n=== RESULTS 72 vs RF: cor ===\n")
print(res_B_cor %>%
        dplyr::select(combo, Higher, pct_x_higher, median_diff,
                      r_Rosenthal, r_rb, p_bonf, sig),
      row.names = FALSE, digits = 4)


cat("\n=== RESULTS 72 vs RF: str_r2 ===\n")
print(res_B_str %>%
        dplyr::select(combo, Higher, pct_x_higher, median_diff,
                      r_Rosenthal, r_rb, p_bonf, sig),
      row.names = FALSE, digits = 4)

# Summary per W_scheme
cat("\n--- Summary Part B per W_scheme ---\n")
for (mname in c("cor", "str_r2")) {
  res_tmp <- if (mname == "cor") res_B_cor else res_B_str
  for (ws in c("robust", "scale")) {
    r_ws <- res_tmp %>% filter(W_scheme == ws)
    n_sword_higher <- sum(r_ws$Higher == "SWORD")
    n_sig          <- sum(r_ws$sig != "n.s.")
    med_r_rb       <- round(median(r_ws$r_rb), 3)
    cat(sprintf("  %s [%s]: SWORD > RF in %d/36 | %d sig | median r_rb = %+.3f\n",
                mname, ws, n_sword_higher, n_sig, med_r_rb))
  }
}

# Global summary
cat("\n--- Global summary Part B (72 configs) ---\n")
for (mname in c("cor", "str_r2")) {
  res_tmp <- if (mname == "cor") res_B_cor else res_B_str
  n_sword_higher <- sum(res_tmp$Higher == "SWORD")
  n_sig          <- sum(res_tmp$sig != "n.s.")
  med_r_rb       <- round(median(res_tmp$r_rb), 3)
  cat(sprintf("  %s: SWORD > RF in %d/72 config | %d sig (p_bonf<0.05) | median r_rb = %+.3f\n",
              mname, n_sword_higher, n_sig, med_r_rb))
}

# =============================================================================
# FORMATTED TABLE - as tab:wilcox_corstr_params in the paper
# Columns: Parameter | Comparison | N | Δ_med(cor) | r_rb(cor) | sig(cor) |
#                                     | Δ_med(str) | r_rb(str) | sig(str)
# =============================================================================

# Readable labels for parameter and comparison
param_labels <- c(
  W_scheme   = "Weights",
  rand_label = "gamma",
  relation   = "Dependency",
  Cost_C     = "C",
  rf_var_frac = "alpha"
)
confront_labels <- c(
  "robust vs scale"   = "robust vs standard",
  "Fixed vs Random"   = "Fixed vs Random",
  "Pearson vs MI"     = "Pearson vs MI",
  "C=0.1 vs C=1"     = "C=0.1 vs C=1",
  "C=0.1 vs C=10"    = "C=0.1 vs C=10",
  "C=1 vs C=10"      = "C=1 vs C=10",
  "30% vs 60%"       = "30% vs 60%",
  "30% vs 100%"      = "30% vs 100%",
  "60% vs 100%"      = "60% vs 100%"
)

# Row order
row_order <- c(
  "Weights|robust vs standard",
  "gamma|Fixed vs Random",
  "Dependency|Pearson vs MI",
  "C|C=0.1 vs C=1",
  "C|C=0.1 vs C=10",
  "C|C=1 vs C=10",
  "alpha|30% vs 60%",
  "alpha|30% vs 100%",
  "alpha|60% vs 100%"
)

# Pivot: cor and str_r2 side by side
tab_wide <- res_A_all %>%
  mutate(
    Parameter  = param_labels[Parametro],
    Comparison = ifelse(Confronto %in% names(confront_labels),
                        confront_labels[Confronto], Confronto),
    row_key    = paste0(Parameter, "|", Comparison)
  ) %>%
  dplyr::select(row_key, Parameter, Comparison, Metrica,
                N_coppie, median_diff, r_rb, sig) %>%
  pivot_wider(
    names_from  = Metrica,
    values_from = c(median_diff, r_rb, sig)
  ) %>%
  mutate(row_key = factor(row_key, levels = row_order)) %>%
  arrange(row_key) %>%
  dplyr::select(
    Parameter, Comparison, N = N_coppie,
    `Δ_med(cor)`  = median_diff_cor,
    `r_rb(cor)`   = r_rb_cor,
    `sig(cor)`    = sig_cor,
    `Δ_med(str)`  = median_diff_str_r2,
    `r_rb(str)`   = r_rb_str_r2,
    `sig(str)`    = sig_str_r2
  )

cat("\n")
cat(paste(rep("=", 90), collapse = ""), "\n")
cat("TABLE: Effect of 5 parameters on Correlation and Strength (R²)\n")
cat("(paper format: tab:wilcox_corstr_params)\n")
cat(paste(rep("=", 90), collapse = ""), "\n")
print(tab_wide, row.names = FALSE, digits = 3)
cat(paste(rep("=", 90), collapse = ""), "\n")

# Also save this table as CSV
write.csv(tab_wide,
          file.path(dir_out, "tab_wilcox_corstr_params_72.csv"),
          row.names = FALSE)

# =============================================================================
# SAVING
# =============================================================================

res_B_all <- bind_rows(res_B_cor, res_B_str)

write.csv(res_A_all,
          file.path(dir_out, "Wilcoxon_corstr_5parametri_72.csv"),
          row.names = FALSE)

write.csv(res_B_all,
          file.path(dir_out, "Wilcoxon_corstr_72vsRF.csv"),
          row.names = FALSE)

cat(sprintf("\nFiles saved in %s:\n  - Wilcoxon_corstr_5parametri_72.csv\n  - Wilcoxon_corstr_72vsRF.csv\n  - tab_wilcox_corstr_params_72.csv\n",
            dir_out))
