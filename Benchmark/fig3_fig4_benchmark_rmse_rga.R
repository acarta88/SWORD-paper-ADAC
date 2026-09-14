###############################################################
# Fig. 3 (fig:bench_res_all, RMSE) + Fig. 4 (fig:bench_res_RGA, RGA)
# BENCHMARK ADAC – Model comparison (with SWORDSCALE)
# Models: SWORDSCALE, RF, ODRF, ROT, EXTree, RRF, aorsfNET, aorsf, SPORF, CF, SVMlin
# Metrics: RMSE, MAE, R2, RGA
# 248 datasets, 50 repetitions
###############################################################

## ============================================================
## Path-robust header — no absolute setwd().
## Works both with `Rscript` and when sourced from RStudio.
##   PROJ    = repo root (reproducibility_code/)
##   DATA    = PROJ/data
##   FIGDIR  = PROJ/figures
##   dir_res = Benchmark/risultati_finali_benchmark_ADAC (next to this script)
## ============================================================
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
DATA    <- file.path(PROJ, "data")
FIGDIR  <- file.path(PROJ, "figures")
CACHE   <- file.path(PROJ, "results_cache")
dir_res <- file.path(SCRIPT_DIR, "risultati_finali_benchmark_ADAC")
if (!dir.exists(FIGDIR)) dir.create(FIGDIR, recursive = TRUE)
if (!interactive()) pdf(NULL)   # avoid Rplots.pdf from on-screen print()

library(dplyr)
library(tidyr)
library(ggplot2)

# ──────────────────────────────────────────────────────────────
# 1-3. LOAD ALL MODELS' RESULTS (SWORDSCALE + 10 competitors) -> df_all
#      Cache-first: use the light cache if present, otherwise rebuild
#      from the per-dataset raw files (11362 files, ~31 MB).
# ──────────────────────────────────────────────────────────────
benchmark_cache_path <- file.path(CACHE, "df_benchmark248_cached.rds")

if (file.exists(benchmark_cache_path)) {

  cat("Loading benchmark results (11 models x 248 datasets) from cache...\n")
  df_all <- readRDS(benchmark_cache_path)
  cat("  Total rows:", nrow(df_all), "\n")
  cat("  Total datasets:", n_distinct(df_all$dataset), "\n")
  cat("  Models:", paste(unique(df_all$model), collapse=", "), "\n\n")

} else {

  cat("Benchmark cache not found - rebuilding from .rds (needs the per-dataset raw files)\n")

  # ---- 1. LOAD SWORDSCALE RESULTS ----
  #      File batch:  risultati_ADAC_SWORDSCALE_<ds>_1_50
  cat("=== Loading SWORDSCALE ===\n")

  swordscale_files_batch <- list.files(dir_res,
                                  pattern = "^risultati_ADAC_SWORDSCALE_.*_[0-9]+_[0-9]+\\.rds$",
                                  full.names = TRUE)
  swordscale_files_rep <- list.files(dir_res,
                                pattern = "^risultati_ADAC_SWORDSCALE_.*_rep[0-9]+\\.rds$",
                                full.names = TRUE)

  # Exclude the aggregated ALL file (if present)
  swordscale_files_batch <- swordscale_files_batch[!grepl("_ALL_", swordscale_files_batch)]

  cat("  File batch:", length(swordscale_files_batch), "\n")
  cat("  Single-rep files:", length(swordscale_files_rep), "\n")

  read_safe <- function(f) {
    tryCatch(readRDS(f), error = function(e) {
      warning("Could not read: ", f); NULL
    })
  }

  list_swordscale_batch <- lapply(swordscale_files_batch, read_safe)
  list_swordscale_rep   <- lapply(swordscale_files_rep,   read_safe)

  df_swordscale_raw <- do.call(rbind, c(list_swordscale_batch, list_swordscale_rep))

  # Add model and spec columns
  df_swordscale_raw$model <- "SWORD"
  df_swordscale_raw$spec  <- "SWORD"

  # Deduplication: keep a single row per (dataset, repetition)
  df_swordscale <- df_swordscale_raw %>%
    arrange(dataset, repetition) %>%
    group_by(dataset, repetition) %>%
    slice_tail(n = 1) %>%
    ungroup()

  cat("  SWORDSCALE rows after dedup:", nrow(df_swordscale), "\n")
  cat("  Dataset SWORDSCALE:", n_distinct(df_swordscale$dataset), "\n")
  cat("  Repetitions per dataset (range):",
      range(table(df_swordscale$dataset)), "\n\n")


  # ---- 2. LOAD OTHER MODELS' RESULTS ----
  cat("=== Loading other models ===\n")

  OTHER_MODELS <- c("RF", "ODRF", "ROT", "EXTree", "RRF",
                    "aorsfNET", "aorsf", "SPORF", "CF", "SVMlin")

  load_model <- function(model_name) {
    pattern <- paste0("^risultati_", model_name, "_ADAC_.*\\.rds$")
    files   <- list.files(dir_res, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) { cat("  No file for:", model_name, "\n"); return(NULL) }
    files <- files[!grepl("_ALL_", files)]
    lst <- lapply(files, read_safe)
    df  <- do.call(rbind, lst)
    if (!is.null(df) && nrow(df) > 0) {
      df <- df[!duplicated(df[, c("dataset", "repetition")]), ]
    }
    if (!"model" %in% names(df)) df$model <- model_name
    if (!"spec"  %in% names(df)) df$spec  <- model_name
    df$model <- model_name
    df$spec  <- model_name
    cat("  ", model_name, ":", nrow(df), "rows,",
        n_distinct(df$dataset), "dataset\n")
    df
  }

  list_other <- lapply(OTHER_MODELS, load_model)
  names(list_other) <- OTHER_MODELS


  # ---- 3. MERGE INTO A SINGLE DATA FRAME ----
  cat("\n=== Merge ===\n")

  CORE_COLS <- c("dataset", "repetition", "RMSE", "MAE", "R2", "RGA",
                 "sd_y_train", "sd_y_test", "ncolumn",
                 "cor_tree", "strength_tree", "time_forest",
                 "model", "spec", "ntree")

  harmonise <- function(df, cols = CORE_COLS) {
    if (is.null(df)) return(NULL)
    for (col in cols) {
      if (!col %in% names(df)) df[[col]] <- NA
    }
    df[, cols, drop = FALSE]
  }

  df_all <- do.call(rbind,
                    c(list(harmonise(df_swordscale)),
                      lapply(list_other, harmonise)))

  cat("  Total rows:", nrow(df_all), "\n")
  cat("  Total datasets:", n_distinct(df_all$dataset), "\n")
  cat("  Models:", paste(unique(df_all$model), collapse=", "), "\n\n")

  if (!dir.exists(CACHE)) dir.create(CACHE, recursive = TRUE)
  saveRDS(df_all, benchmark_cache_path)
  cat("Cache saved to", benchmark_cache_path, "\n\n")

}


# ──────────────────────────────────────────────────────────────
# 4. RMSE NORMALIZATION (NRMSE = RMSE / sd_y_test)
# ──────────────────────────────────────────────────────────────
df_all <- df_all %>%
  mutate(
    NRMSE = RMSE / sd_y_test,
    NMAE  = MAE  / sd_y_test
  )



# ──────────────────────────────────────────────────────────────
# 5. IMPROVEMENT vs RF (computed at the single-repetition level)
# ──────────────────────────────────────────────────────────────
cat("=== Computing improvement vs RF (per repetition) ===\n")

df_unito_248 <- df_all %>%
  group_by(dataset, repetition) %>%
  mutate(
    RMSE_RF          = if (any(model == "RF")) RMSE[model == "RF"][1] else NA_real_,
    RGA_RF           = if (any(model == "RF")) RGA [model == "RF"][1] else NA_real_,
    improvement_RMSE = (RMSE_RF - RMSE) / RMSE_RF,
    improvement_RGA  = (RGA  - RGA_RF) 
  ) %>%
  ungroup()

df_unito <- df_unito_248 %>%
  group_by(dataset, repetition) %>%
  mutate(
    RMSE_RF = if_else(
      any(spec == "RF"),
      RMSE[spec == "RF"][1],
      NA_real_
    ),
    improvement_RF = if_else(
      !is.na(RMSE_RF),
      (RMSE_RF - RMSE) / RMSE_RF,
      NA_real_
    )
  ) %>%
  ungroup()


# ──────────────────────────────────────────────────────────────
# 4.3 Load dataset metadata
# ──────────────────────────────────────────────────────────────
load(file.path(DATA, "Info_df_ADAC.rda"))   # provides Info_248_df

df_unito_enriched <- df_unito %>%
  inner_join(Info_248_df, by = "dataset")


############################################################
# 5. FIGURE – IMPROVEMENT vs RF ACROSS MODELS (SCATTER)
############################################################

ordine_spec <- c(
  "SWORD", "ODRF", "EXTree", "RF", "RRF",
  "aorsfNET", "aorsf", "SPORF", "CF", "ROT", "SVMlin"
)

df_unito_enriched$spec <- factor(
  df_unito_enriched$spec,
  levels = ordine_spec
)

df_media_per_dataset <- df_unito_enriched %>%
  group_by(dataset, spec) %>%
  summarise(
    NRMSE_RF = mean(improvement_RF, na.rm = TRUE),
    .groups  = "drop"
  )

spec_levels <- c(
  "SWORD", "ODRF", "EXTree", "RF", "RRF",
  "aorsfNET", "aorsf", "SPORF", "CF", "ROT", "SVMlin"
)

df_media_per_dataset <- df_media_per_dataset %>%
  dplyr::mutate(spec = factor(spec, levels = spec_levels))

df_media_per_dataset <- df_media_per_dataset[df_media_per_dataset$spec != "RF", ]

appiatimento <- 0.5

dataset_scatter_all <- ggplot(df_media_per_dataset, aes(x = spec, y = NRMSE_RF)) +
  geom_jitter(
    data = df_media_per_dataset %>%
      mutate(
        y_jitter = ifelse(
          NRMSE_RF >  appiatimento,  appiatimento,
          ifelse(NRMSE_RF < -appiatimento, -appiatimento, NRMSE_RF)
        )
      ),
    aes(y = y_jitter),
    width  = 0.2,
    height = 0,
    alpha  = 0.5,
    size   = 1.5,
    color  = "darkgray"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  stat_summary(fun = mean,   geom = "point", shape = 4, size = 2, stroke = 1, color = "black") +
  labs(title = NULL, x = "Model", y = "All") +
  coord_cartesian(ylim = c(-appiatimento - 0.1, appiatimento + 0.1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dataset_scatter_all


############################################################
# 6. RANKING AND HEATMAP (RMSE)
############################################################

tab_to_perc <- function(df, var_spec, var_rank) {
  tab <- table(df[[var_spec]], df[[var_rank]])
  as.data.frame(tab) %>%
    group_by(Var1) %>%
    mutate(perc = Freq / sum(Freq) * 100) %>%
    ungroup() %>%
    rename(spec = Var1, rank = Var2, count = Freq)
}

df_rank_rmse <- df_unito_enriched %>%
  group_by(dataset, spec) %>%
  summarise(
    mean_RMSE            = mean(RMSE,           na.rm = TRUE),
    mean_R2              = mean(R2,             na.rm = TRUE),
    mean_improvement_RF  = mean(improvement_RF, na.rm = TRUE),
    median_improvement_RF = median(improvement_RF, na.rm = TRUE),
    median_RMSE          = median(RMSE,         na.rm = TRUE),
    sd_RMSE              = sd(RMSE,             na.rm = TRUE),
    .groups              = "drop"
  ) %>%
  group_by(dataset) %>%
  mutate(
    RMSE_rank           = rank(mean_RMSE,          ties.method = "first"),
    improvement_RF_rank = rank(-mean_improvement_RF, ties.method = "first"),
    R2_rank             = rank(-mean_R2,            ties.method = "first")
  ) %>%
  arrange(dataset, RMSE_rank)

df_rank_rmse_enriched <- df_rank_rmse %>%
  inner_join(Info_248_df, by = "dataset")

df_all_perc <- tab_to_perc(df_rank_rmse_enriched, "spec", "RMSE_rank")
df_all_perc <- df_all_perc %>%
  mutate(spec = factor(spec, levels = rev(spec_levels)))

heatmap_all <- ggplot(df_all_perc, aes(x = rank, y = spec, fill = perc)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "black", name = "% Freq.", limits = c(0, 60)) +
  labs(title = NULL, x = "Rank", y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

heatmap_all


############################################################
# 7. CATEGORICAL vs NUMERICAL – SCATTER + HEATMAP (RMSE)
############################################################

df_media_per_dataset_enriched <- df_media_per_dataset %>%
  left_join(Info_248_df, by = "dataset")
df_media_per_dataset_enriched <- df_media_per_dataset_enriched[
  df_media_per_dataset_enriched$spec != "RF", ]

# Categorical
df_media_per_dataset_enriched_cat <- df_media_per_dataset_enriched[
  df_media_per_dataset_enriched$n_categorical_features != 0 |
    df_media_per_dataset_enriched$n_binary_features != 0, ]

df_rank_rmse_enriched_cat <- df_rank_rmse_enriched[
  df_rank_rmse_enriched$n_categorical_features != 0 |
    df_rank_rmse_enriched$n_binary_features != 0, ]

dataset_scatter_cat <- ggplot(df_media_per_dataset_enriched_cat, aes(x = spec, y = NRMSE_RF)) +
  geom_jitter(
    data = df_media_per_dataset_enriched_cat %>%
      mutate(y_jitter = ifelse(NRMSE_RF > appiatimento, appiatimento,
                               ifelse(NRMSE_RF < -appiatimento, -appiatimento, NRMSE_RF))),
    aes(y = y_jitter), width = 0.2, height = 0, alpha = 0.5, size = 1.5, color = "darkgray"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  stat_summary(fun = mean,   geom = "point", shape = 4, size = 2, stroke = 1, color = "black") +
  labs(title = NULL, x = NULL, y = "Categorical") +
  coord_cartesian(ylim = c(-appiatimento - 0.1, appiatimento + 0.1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

df_cat_perc <- tab_to_perc(df_rank_rmse_enriched_cat, "spec", "RMSE_rank")
heatmap_cat <- ggplot(df_cat_perc, aes(x = rank, y = spec, fill = perc)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "black", name = "% Freq.", limits = c(0, 60)) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(limits = rev(spec_levels))

# Numerical
df_media_per_dataset_enriched_num <- df_media_per_dataset_enriched[
  df_media_per_dataset_enriched$n_categorical_features == 0 &
    df_media_per_dataset_enriched$n_binary_features == 0, ]

df_rank_rmse_enriched_num <- df_rank_rmse_enriched[
  df_rank_rmse_enriched$n_categorical_features == 0 &
    df_rank_rmse_enriched$n_binary_features == 0, ]

dataset_scatter_num <- ggplot(df_media_per_dataset_enriched_num, aes(x = spec, y = NRMSE_RF)) +
  geom_jitter(
    data = df_media_per_dataset_enriched_num %>%
      mutate(y_jitter = ifelse(NRMSE_RF > appiatimento, appiatimento,
                               ifelse(NRMSE_RF < -appiatimento, -appiatimento, NRMSE_RF))),
    aes(y = y_jitter), width = 0.2, height = 0, alpha = 0.5, size = 1.5, color = "darkgray"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  stat_summary(fun = mean,   geom = "point", shape = 4, size = 2, stroke = 1, color = "black") +
  labs(title = NULL, x = NULL, y = "Numeric") +
  coord_cartesian(ylim = c(-appiatimento - 0.1, appiatimento + 0.1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

df_num_perc <- tab_to_perc(df_rank_rmse_enriched_num, "spec", "RMSE_rank")
heatmap_num <- ggplot(df_num_perc, aes(x = rank, y = spec, fill = perc)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "black", name = "% Freq.", limits = c(0, 60)) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(limits = rev(spec_levels))


############################################################
# 8. COMBINED FIGURE – SCATTER + HEATMAP (2×3 PANEL) – RMSE
############################################################

library(patchwork)
library(cowplot)

col_a <- (dataset_scatter_num / dataset_scatter_cat / dataset_scatter_all) +
  plot_layout(heights = c(1, 1, 1)) &
  theme(legend.position = "none")

col_b <- (heatmap_num / heatmap_cat / heatmap_all) +
  plot_layout(guides = "collect", heights = c(1, 1, 1)) &
  theme(legend.position = "right")

combined_2x3 <- cowplot::plot_grid(
  col_a, col_b,
  labels     = c("a", "b"),
  label_size = 14,
  ncol       = 2,
  align      = "h",
  rel_widths = c(1, 1.08)
)

combined_2x3


##############################  RGA
##############################  RGA
##############################  RGA

MODEL_ORDER <- c("SWORD", "RF", "ODRF", "EXTree", "RRF",
                 "aorsfNET", "aorsf", "SPORF", "CF", "ROT", "SVMlin")

summary_by_ds <- df_unito_enriched %>%
  group_by(model, dataset) %>%
  summarise(
    n_rep          = n(),
    mean_RMSE      = mean(RMSE,             na.rm = TRUE),
    mean_NRMSE     = mean(NRMSE,            na.rm = TRUE),
    mean_MAE       = mean(MAE,              na.rm = TRUE),
    mean_R2        = mean(R2,               na.rm = TRUE),
    mean_RGA       = mean(RGA,              na.rm = TRUE),
    mean_time      = mean(time_forest,      na.rm = TRUE),
    mean_impRMSE   = mean(improvement_RMSE, na.rm = TRUE),
    mean_impRGA    = mean(improvement_RGA,  na.rm = TRUE),
    .groups = "drop"
  )

summary_models <- summary_by_ds %>%
  group_by(model) %>%
  summarise(
    n_dataset        = n(),
    avg_NRMSE        = mean(mean_NRMSE,   na.rm = TRUE),
    med_NRMSE        = median(mean_NRMSE, na.rm = TRUE),
    avg_R2           = mean(mean_R2,      na.rm = TRUE),
    med_R2           = median(mean_R2,    na.rm = TRUE),
    avg_RGA          = mean(mean_RGA,     na.rm = TRUE),
    med_RGA          = median(mean_RGA,   na.rm = TRUE),
    avg_time         = mean(mean_time,    na.rm = TRUE),
    avg_impRMSE      = mean(mean_impRMSE,   na.rm = TRUE),
    med_impRMSE      = median(mean_impRMSE, na.rm = TRUE),
    pct_beat_RF_RMSE = mean(mean_impRMSE > 0, na.rm = TRUE) * 100,
    avg_impRGA       = mean(mean_impRGA,   na.rm = TRUE),
    med_impRGA       = median(mean_impRGA, na.rm = TRUE),
    pct_beat_RF_RGA  = mean(mean_impRGA > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(model = factor(model, levels = MODEL_ORDER)) %>%
  arrange(model)

cat("\n=== FINAL SUMMARY (mean and median across datasets) ===\n")
print(as.data.frame(summary_models), digits = 4)


# ──────────────────────────────────────────────────────────────
# TABLE: mean, median, SD of RMSE and RGA improvement per model
# ──────────────────────────────────────────────────────────────

tabella_improvement <- summary_by_ds %>%
  group_by(model) %>%
  summarise(
    n_dataset          = n(),
    # Improvement RMSE vs RF
    mean_imp_RMSE      = mean(mean_impRMSE,    na.rm = TRUE),
    median_imp_RMSE    = median(mean_impRMSE,  na.rm = TRUE),
    sd_imp_RMSE        = sd(mean_impRMSE,      na.rm = TRUE),
    # Improvement RGA vs RF
    mean_imp_RGA       = mean(mean_impRGA,     na.rm = TRUE),
    median_imp_RGA     = median(mean_impRGA,   na.rm = TRUE),
    sd_imp_RGA         = sd(mean_impRGA,       na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(model = factor(model, levels = MODEL_ORDER)) %>%
  arrange(model)

cat("\n=== TABLE: Mean / Median / SD of RMSE and RGA Improvement per model ===\n")
cat("(Improvement RMSE = (RMSE_RF - RMSE) / RMSE_RF, positive = better than RF)\n")
cat("(Improvement RGA  = RGA - RGA_RF,                positive = better than RF)\n\n")
print(as.data.frame(tabella_improvement), digits = 4)

df_media_per_dataset_RGA <- df_unito_enriched %>%
  group_by(dataset, spec) %>%
  summarise(
    RGA_PERC_RF = mean(improvement_RGA, na.rm = TRUE),
    .groups  = "drop"
  )

df_media_per_dataset_RGA <- df_media_per_dataset_RGA %>%
  dplyr::mutate(spec = factor(spec, levels = spec_levels))
df_media_per_dataset_RGA <- df_media_per_dataset_RGA[df_media_per_dataset_RGA$spec != "RF", ]

appiatimento_rga <- 0.05

dataset_scatter_RGA_all <- ggplot(df_media_per_dataset_RGA, aes(x = spec, y = RGA_PERC_RF)) +
  geom_jitter(
    data = df_media_per_dataset_RGA %>%
      mutate(y_jitter = ifelse(RGA_PERC_RF > appiatimento_rga, appiatimento_rga,
                               ifelse(RGA_PERC_RF < -appiatimento_rga, -appiatimento_rga, RGA_PERC_RF))),
    aes(y = y_jitter), width = 0.2, height = 0, alpha = 0.5, size = 1.5, color = "darkgray"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  stat_summary(fun = mean,   geom = "point", shape = 4, size = 2, stroke = 1, color = "black") +
  labs(title = NULL, x = "Model", y = "All") +
  coord_cartesian(ylim = c(-appiatimento_rga - 0.01, appiatimento_rga + 0.01)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dataset_scatter_RGA_all

df_rank_RGA <- df_unito_enriched %>%
  group_by(dataset, spec) %>%
  summarise(
    mean_RGA                 = mean(RGA,             na.rm = TRUE),
    mean_R2                  = mean(R2,              na.rm = TRUE),
    mean_improvement_RGA_RF  = mean(improvement_RGA, na.rm = TRUE),
    median_improvement_RGA_RF = median(improvement_RGA, na.rm = TRUE),
    median_RGA               = median(RGA,           na.rm = TRUE),
    sd_RGA                   = sd(RGA,               na.rm = TRUE),
    .groups                  = "drop"
  ) %>%
  group_by(dataset) %>%
  mutate(
    RGA_rank               = rank(-mean_RGA,              ties.method = "first"),
    improvement_RGA_RF_rank = rank(-mean_improvement_RGA_RF, ties.method = "first"),
    R2_rank                = rank(-mean_R2,               ties.method = "first")
  ) %>%
  arrange(dataset, RGA_rank)

df_rank_RGA_enriched <- df_rank_RGA %>%
  inner_join(Info_248_df, by = "dataset")

df_all_rga_perc <- tab_to_perc(df_rank_RGA_enriched, "spec", "RGA_rank")
df_all_rga_perc <- df_all_rga_perc %>%
  mutate(spec = factor(spec, levels = rev(spec_levels)))

heatmap_RGA_all <- ggplot(df_all_rga_perc, aes(x = rank, y = spec, fill = perc)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "black", name = "% Freq.", limits = c(0, 60)) +
  labs(title = NULL, x = "Rank", y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

heatmap_RGA_all

df_media_per_dataset_RGA_enriched <- df_media_per_dataset_RGA %>%
  left_join(Info_248_df, by = "dataset")
df_media_per_dataset_RGA_enriched <- df_media_per_dataset_RGA_enriched[
  df_media_per_dataset_RGA_enriched$spec != "RF", ]

# Categorical – RGA
df_media_per_dataset_RGA_enriched_cat <- df_media_per_dataset_RGA_enriched[
  df_media_per_dataset_RGA_enriched$n_categorical_features != 0 |
    df_media_per_dataset_RGA_enriched$n_binary_features != 0, ]
df_rank_RGA_enriched_cat <- df_rank_RGA_enriched[
  df_rank_RGA_enriched$n_categorical_features != 0 |
    df_rank_RGA_enriched$n_binary_features != 0, ]

dataset_scatter_RGA_cat <- ggplot(df_media_per_dataset_RGA_enriched_cat, aes(x = spec, y = RGA_PERC_RF)) +
  geom_jitter(
    data = df_media_per_dataset_RGA_enriched_cat %>%
      mutate(y_jitter = ifelse(RGA_PERC_RF > appiatimento_rga, appiatimento_rga,
                               ifelse(RGA_PERC_RF < -appiatimento_rga, -appiatimento_rga, RGA_PERC_RF))),
    aes(y = y_jitter), width = 0.2, height = 0, alpha = 0.5, size = 1.5, color = "darkgray"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  stat_summary(fun = mean,   geom = "point", shape = 4, size = 2, stroke = 1, color = "black") +
  labs(title = NULL, x = NULL, y = "Categorical") +
  coord_cartesian(ylim = c(-appiatimento_rga - 0.01, appiatimento_rga + 0.01)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

df_cat_rga_perc <- tab_to_perc(df_rank_RGA_enriched_cat, "spec", "RGA_rank")
heatmap_RGA_cat <- ggplot(df_cat_rga_perc, aes(x = rank, y = spec, fill = perc)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "black", name = "% Freq.", limits = c(0, 60)) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(limits = rev(spec_levels))

# Numerical – RGA
df_media_per_dataset_RGA_enriched_num <- df_media_per_dataset_RGA_enriched[
  df_media_per_dataset_RGA_enriched$n_categorical_features == 0 &
    df_media_per_dataset_RGA_enriched$n_binary_features == 0, ]
df_rank_RGA_enriched_num <- df_rank_RGA_enriched[
  df_rank_RGA_enriched$n_categorical_features == 0 &
    df_rank_RGA_enriched$n_binary_features == 0, ]

dataset_scatter_RGA_num <- ggplot(df_media_per_dataset_RGA_enriched_num, aes(x = spec, y = RGA_PERC_RF)) +
  geom_jitter(
    data = df_media_per_dataset_RGA_enriched_num %>%
      mutate(y_jitter = ifelse(RGA_PERC_RF > appiatimento_rga, appiatimento_rga,
                               ifelse(RGA_PERC_RF < -appiatimento_rga, -appiatimento_rga, RGA_PERC_RF))),
    aes(y = y_jitter), width = 0.2, height = 0, alpha = 0.5, size = 1.5, color = "darkgray"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  stat_summary(fun = mean,   geom = "point", shape = 4, size = 2, stroke = 1, color = "black") +
  labs(title = NULL, x = NULL, y = "Numeric") +
  coord_cartesian(ylim = c(-appiatimento_rga - 0.01, appiatimento_rga + 0.01)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

df_num_rga_perc <- tab_to_perc(df_rank_RGA_enriched_num, "spec", "RGA_rank")
heatmap_RGA_num <- ggplot(df_num_rga_perc, aes(x = rank, y = spec, fill = perc)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "black", name = "% Freq.", limits = c(0, 60)) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(limits = rev(spec_levels))


############################################################
# 8. COMBINED FIGURE – SCATTER + HEATMAP (2×3 PANEL) – RGA
############################################################

col_a_RGA <- (dataset_scatter_RGA_num / dataset_scatter_RGA_cat / dataset_scatter_RGA_all) +
  plot_layout(heights = c(1, 1, 1)) &
  theme(legend.position = "none")

col_b_RGA <- (heatmap_RGA_num / heatmap_RGA_cat / heatmap_RGA_all) +
  plot_layout(guides = "collect", heights = c(1, 1, 1)) &
  theme(legend.position = "right")

combined_2x3_RGA <- cowplot::plot_grid(
  col_a_RGA, col_b_RGA,
  labels     = c("a", "b"),
  label_size = 14,
  ncol       = 2,
  align      = "h",
  rel_widths = c(1, 1.08)
)

combined_2x3_RGA
combined_2x3


############################################################
# IMPROVEMENT TABLE – NUMERIC VS CATEGORICAL
############################################################

tabella_improvement_tipo <- summary_by_ds %>%
  left_join(
    Info_248_df %>%
      select(dataset, n_categorical_features, n_binary_features) %>%
      mutate(
        data_type = ifelse(
          n_categorical_features != 0 | n_binary_features != 0,
          "Categorical",
          "Numeric"
        )
      ) %>%
      select(dataset, data_type),
    by = "dataset"
  ) %>%
  group_by(data_type, model) %>%
  summarise(
    n_dataset          = n(),
    # Improvement RMSE vs RF
    mean_imp_RMSE      = mean(mean_impRMSE,   na.rm = TRUE),
    median_imp_RMSE    = median(mean_impRMSE, na.rm = TRUE),
    sd_imp_RMSE        = sd(mean_impRMSE,     na.rm = TRUE),
    # Improvement RGA vs RF
    mean_imp_RGA       = mean(mean_impRGA,    na.rm = TRUE),
    median_imp_RGA     = median(mean_impRGA,  na.rm = TRUE),
    sd_imp_RGA         = sd(mean_impRGA,      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    data_type = factor(data_type, levels = c("Numeric", "Categorical")),
    model     = factor(model, levels = MODEL_ORDER)
  ) %>%
  arrange(data_type, model)

cat("\n=== TABLE: Mean / Median / SD of RMSE and RGA Improvement per model and dataset type ===\n")
cat("(Improvement RMSE = (RMSE_RF - RMSE) / RMSE_RF, positive = better than RF)\n")
cat("(Improvement RGA  = RGA - RGA_RF,                positive = better than RF)\n\n")
print(as.data.frame(tabella_improvement_tipo), digits = 4)
print(as.data.frame(tabella_improvement), digits = 4)


############################################################
# TABLE: Rank 1, Rank 2, Beats RF - per metric and dataset type
############################################################

# Helper: adds data_type column to a df already joined with Info_248_df
add_data_type <- function(df) {
  df %>% mutate(
    data_type = ifelse(
      n_categorical_features != 0 | n_binary_features != 0,
      "Categorical", "Numeric"
    )
  )
}

# ── RMSE ──────────────────────────────────────────────────────────────────────

rmse_all <- df_rank_rmse_enriched %>%
  group_by(spec) %>%
  summarise(
    n_dataset      = n(),
    n_rank1_RMSE   = sum(RMSE_rank == 1, na.rm = TRUE),
    n_rank2_RMSE   = sum(RMSE_rank == 2, na.rm = TRUE),
    n_beat_RF_RMSE = sum(mean_improvement_RF > 0, na.rm = TRUE),
    data_type      = "All",
    .groups = "drop"
  )

rmse_tipo <- add_data_type(df_rank_rmse_enriched) %>%
  group_by(spec, data_type) %>%
  summarise(
    n_dataset      = n(),
    n_rank1_RMSE   = sum(RMSE_rank == 1, na.rm = TRUE),
    n_rank2_RMSE   = sum(RMSE_rank == 2, na.rm = TRUE),
    n_beat_RF_RMSE = sum(mean_improvement_RF > 0, na.rm = TRUE),
    .groups = "drop"
  )

tab_rmse <- bind_rows(rmse_all, rmse_tipo)

# ── RGA ───────────────────────────────────────────────────────────────────────

rga_all <- df_rank_RGA_enriched %>%
  group_by(spec) %>%
  summarise(
    n_dataset     = n(),
    n_rank1_RGA   = sum(RGA_rank == 1, na.rm = TRUE),
    n_rank2_RGA   = sum(RGA_rank == 2, na.rm = TRUE),
    n_beat_RF_RGA = sum(mean_improvement_RGA_RF > 0, na.rm = TRUE),
    data_type     = "All",
    .groups = "drop"
  )

rga_tipo <- add_data_type(df_rank_RGA_enriched) %>%
  group_by(spec, data_type) %>%
  summarise(
    n_dataset     = n(),
    n_rank1_RGA   = sum(RGA_rank == 1, na.rm = TRUE),
    n_rank2_RGA   = sum(RGA_rank == 2, na.rm = TRUE),
    n_beat_RF_RGA = sum(mean_improvement_RGA_RF > 0, na.rm = TRUE),
    .groups = "drop"
  )

tab_rga <- bind_rows(rga_all, rga_tipo)

# -- Merge ────────────────────────────────────────────────────────────────────

tab_rank_summary <- tab_rmse %>%
  full_join(tab_rga %>% select(-n_dataset), by = c("spec", "data_type")) %>%
  mutate(
    data_type = factor(data_type, levels = c("All", "Numeric", "Categorical")),
    spec      = factor(spec, levels = spec_levels)
  ) %>%
  arrange(data_type, spec) %>%
  rename(model = spec)

tab_rank_pct <- tab_rank_summary %>%
  mutate(
    pct_rank1_RMSE   = round(n_rank1_RMSE   / n_dataset * 100, 1),
    pct_rank2_RMSE   = round(n_rank2_RMSE   / n_dataset * 100, 1),
    pct_beat_RF_RMSE = round(n_beat_RF_RMSE / n_dataset * 100, 1),
    pct_rank1_RGA    = round(n_rank1_RGA    / n_dataset * 100, 1),
    pct_rank2_RGA    = round(n_rank2_RGA    / n_dataset * 100, 1),
    pct_beat_RF_RGA  = round(n_beat_RF_RGA  / n_dataset * 100, 1),
    # Labels "n (x%)" for readability
    rank1_RMSE   = paste0(n_rank1_RMSE,   " (", pct_rank1_RMSE,   "%)"),
    rank2_RMSE   = paste0(n_rank2_RMSE,   " (", pct_rank2_RMSE,   "%)"),
    beat_RF_RMSE = paste0(n_beat_RF_RMSE, " (", pct_beat_RF_RMSE, "%)"),
    rank1_RGA    = paste0(n_rank1_RGA,    " (", pct_rank1_RGA,    "%)"),
    rank2_RGA    = paste0(n_rank2_RGA,    " (", pct_rank2_RGA,    "%)"),
    beat_RF_RGA  = paste0(n_beat_RF_RGA,  " (", pct_beat_RF_RGA,  "%)")
  )

cat("\n=== TABLE: Rank 1, Rank 2, Beats RF - n (%) per metric and dataset type ===\n")
for (tipo in c("All", "Numeric", "Categorical")) {
  cat(sprintf("\n--- %s (n = %d) ---\n", tipo,
              tab_rank_pct$n_dataset[tab_rank_pct$data_type == tipo][1]))
  sub <- tab_rank_pct %>%
    filter(data_type == tipo) %>%
    select(model,
           rank1_RMSE, rank2_RMSE, beat_RF_RMSE,
           rank1_RGA,  rank2_RGA,  beat_RF_RGA)
  print(as.data.frame(sub), row.names = FALSE)
}


############################################################
# SAVING THE PAPER'S FIGURES
#   fig:bench_res_all  -> bench_2x3_panels_rmse_svmlin.pdf
#   fig:bench_res_RGA  -> bench_2x3_panels_RGA_svmlin.pdf
############################################################

ggsave(file.path(FIGDIR, "bench_2x3_panels_rmse_svmlin.pdf"),
       combined_2x3,     width = 8, height = 9.3)
ggsave(file.path(FIGDIR, "bench_2x3_panels_RGA_svmlin.pdf"),
       combined_2x3_RGA, width = 8, height = 9.3)

cat("\nFigures saved in:", FIGDIR, "\n")
cat("  - bench_2x3_panels_rmse_svmlin.pdf\n")
cat("  - bench_2x3_panels_RGA_svmlin.pdf\n")

