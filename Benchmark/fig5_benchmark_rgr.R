# =============================================================================
# Fig. 5 (fig:bench_rgr) — RGR / delta-RGR figure for the paper -> fig_RGR_delta_2x2.pdf
# Figure stage only: reads the RGR_<model>.rds files already computed in
# Benchmark/output_RGR (no recomputation from raw predictions).
# =============================================================================

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
DATA    <- file.path(PROJ, "data")
FIGDIR  <- file.path(PROJ, "figures")
TABDIR  <- file.path(PROJ, "tables")
OUT_RGR <- file.path(SCRIPT_DIR, "output_RGR")
if (!dir.exists(FIGDIR)) dir.create(FIGDIR, recursive = TRUE)
if (!dir.exists(TABDIR)) dir.create(TABDIR, recursive = TRUE)
if (!interactive()) pdf(NULL)

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(cowplot)

# =============================================================================
# 1. LOAD DATA
# =============================================================================


out_dir <- OUT_RGR

# Load all RGR_<model>.rds files (exclude _state.rds)
files <- list.files(out_dir, pattern = "^RGR_.*\\.rds$", full.names = TRUE)
files <- files[!grepl("_state", basename(files))]
cat("Files found:", length(files), "\n")
cat(paste(" -", basename(files)), sep = "\n")

# Read and combine
df_long <- bind_rows(lapply(files, readRDS))
cat("\nTotal rows:", nrow(df_long),
    "| Models:", n_distinct(df_long$model),
    "| Datasets:", n_distinct(df_long$dataset),
    "| Rep:", n_distinct(df_long$rep), "\n")

# Rename SWORDSCALE -> SWORD for consistency
df_long$model[df_long$model == "SWORDSCALE"] <- "SWORD"
df_long <- df_long %>% distinct()
# Average per dataset x model (mean over repetitions)
df_rgr <- df_long %>%
  group_by(model, dataset) %>%
  summarise(mean_RGR = mean(RGR, na.rm = TRUE),
            median_RGR = mean(RGR, na.rm = TRUE),
            sd_RGR   = sd(RGR,   na.rm = TRUE),
            n_rep    = sum(!is.na(RGR)),
            .groups  = "drop")


load(file.path(DATA, "Info_df_ADAC.rda"))            # -> Info_248_df
load(file.path(DATA, "benchmark_248_dataset_list_ADAC.rda"))
# Sanitize dataset names for the join (as in the original file)
Info_248_df_safe <- Info_248_df %>%
  mutate(dataset = gsub("[^A-Za-z0-9_]", "_", dataset))

# =============================================================================
# 2. IMPROVEMENT RGR vs RF  (mean_RGR - RGR_RF per dataset)
# =============================================================================

df_rgr <- df_rgr %>%
  group_by(dataset) %>%
  mutate(
    RGR_RF          = mean_RGR[model == "RF"][1],
    improvement_RGR = mean_RGR - RGR_RF
  ) %>%
  ungroup()

# =============================================================================
# 3. DATASET METADATA (Numeric / Categorical)
# =============================================================================

df_rgr <- df_rgr %>%
  left_join(
    Info_248_df_safe %>%
      select(dataset, n_categorical_features, n_binary_features,n_continuous_features, n_features) %>%
      mutate(
        data_type = ifelse(
          n_categorical_features != 0 | n_binary_features != 0,
          "Categorical", "Numeric"
        )
      ),
    by = "dataset"
  )

dataset_da_togliere <- names(
  Filter(function(df) {
    
    x <- df[, setdiff(names(df), "target"), drop = FALSE]
    
    if (any(sapply(x, is.factor) | sapply(x, is.character))) {
      x <- fastDummies::dummy_cols(
        x,
        remove_selected_columns = TRUE,
        remove_first_dummy = TRUE
      )
    }
    
    all(sapply(x, function(col) length(unique(col[!is.na(col)])) <= 10))
    
  }, benchmark_248_dataset_list)
)

dataset_da_togliere


# remove all datasets whose variables all have fewer than 11 levels
df_rgr <- df_rgr %>%
  filter(!dataset %in% dataset_da_togliere)

cat(sprintf("Rows df_rgr: %d | Datasets: %d | Models: %s\n",
            nrow(df_rgr),
            n_distinct(df_rgr$dataset),
            paste(unique(df_rgr$model), collapse = ", ")))
df_rgr=na.omit(df_rgr)
df_rgr$model[df_rgr$model == "SWORDSCALE"] <- "SWORD"

# =============================================================================
# 4. SUMMARY TABLE
# =============================================================================

MODEL_ORDER <- c("SWORD", "RF", "ODRF", "EXTree", "RRF",
                 "aorsfNET", "aorsf", "SPORF", "CF", "ROT", "SVMlin")

summarise_rgr <- function(df_sub) {
  df_sub %>%
    mutate(model = factor(model, levels = MODEL_ORDER)) %>%
    group_by(model) %>%
    summarise(
      n_dataset      = n(),
      mean_RGR2       = mean(mean_RGR,          na.rm = TRUE),
      median_RGR     = median(mean_RGR,        na.rm = TRUE),
      sd_RGR         = (var(mean_RGR,            na.rm = TRUE))^0.5,
      mean_imp_RGR   = mean(improvement_RGR,   na.rm = TRUE),
      median_imp_RGR = median(improvement_RGR, na.rm = TRUE),
      sd_imp_RGR     = sd(improvement_RGR,     na.rm = TRUE),
      pct_above_0    = mean(improvement_RGR > 0, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    arrange(model)
}

tab_all <- summarise_rgr(df_rgr) %>% mutate(data_type = "All")
tab_num <- summarise_rgr(df_rgr %>% filter(data_type == "Numeric")) %>%
  mutate(data_type = "Numeric")
tab_cat <- summarise_rgr(df_rgr %>% filter(data_type == "Categorical")) %>%
  mutate(data_type = "Categorical")

tabella_RGR <- bind_rows(tab_all, tab_num, tab_cat) %>%
  mutate(data_type = factor(data_type, levels = c("All", "Numeric", "Categorical"))) %>%
  arrange(data_type, model)

cat("\n=== RGR TABLE — by model and dataset type ===\n")
cat("(improvement_RGR = mean_RGR - RGR_RF; >0 = more robust than RF)\n\n")
for (tipo in c("All", "Numeric", "Categorical")) {
  cat(sprintf("\n--- %s ---\n", tipo))
  print(as.data.frame(tabella_RGR %>% filter(data_type == tipo) %>%
                        select(-data_type)), digits = 4, row.names = FALSE)
}

# =============================================================================
# 5. COMMON HELPERS
# =============================================================================

tab_to_perc <- function(df, var_spec, var_rank) {
  tab <- table(df[[var_spec]], df[[var_rank]])
  as.data.frame(tab) %>%
    group_by(Var1) %>%
    mutate(perc = Freq / sum(Freq) * 100) %>%
    ungroup() %>%
    rename(spec = Var1, rank = Var2, count = Freq)
}

# Model order (with RF for absolute RGR)
spec_levels_rgr <- c("SWORD", "ODRF", "EXTree", "RF", "RRF",
                     "aorsfNET", "aorsf", "SPORF", "CF", "ROT", "SVMlin")

# Model order (without RF for delta-RGR)
spec_levels_imp <- c("SWORD", "ODRF", "EXTree", "RRF",
                     "aorsfNET", "aorsf", "SPORF", "CF", "ROT", "SVMlin")

appiatimento_rgr <- 0.15   # clip scatter RGR around 0.5 +/- 0.15 -> [0.35, 0.65]... no, better lateral clip
# For absolute RGR: axis from 0 to 1, lower clip at 0 (no strong capping needed)
# For delta-RGR: clip +/- appiatimento_imp
appiatimento_imp <- 0.15

# =============================================================================
# 6. FIGURE A — ABSOLUTE RGR (2x3 panel: scatter | heatmap)
# =============================================================================

# ── 6a. Rank by RGR (rank 1 = highest value = most robust) ──────────────────

df_rank_rgr <- df_rgr %>%
  mutate(model = factor(model, levels = spec_levels_rgr)) %>%
  group_by(dataset) %>%
  mutate(
    RGR_rank = rank(-mean_RGR, ties.method = "first")
  ) %>%
  ungroup()

df_rank_rgr_enriched <- df_rank_rgr %>%
  left_join(
    Info_248_df_safe %>%
      select(dataset, n_categorical_features, n_binary_features) %>%
      mutate(
        data_type_ = ifelse(
          n_categorical_features != 0 | n_binary_features != 0,
          "Categorical", "Numeric"
        )
      ),
    by = "dataset"
  )

# ── 6b. Scatter RGR by stratum ────────────────────────────────────────────────

make_scatter_rgr <- function(df_sub, y_label, show_xlab = FALSE) {
  df_sub <- df_sub %>%
    filter(model %in% spec_levels_rgr) %>%
    mutate(model = factor(model, levels = spec_levels_rgr))
  
  ggplot(df_sub, aes(x = model, y = mean_RGR)) +
    geom_jitter(
      data = df_sub %>%
        mutate(y_jitter = pmin(pmax(mean_RGR, 0), 1)),
      aes(y = y_jitter),
      width = 0.2, height = 0,
      alpha = 0.5, size = 1.5, color = "darkgray"
    ) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
    stat_summary(fun = mean, geom = "point",
                 shape = 4, size = 2, stroke = 1, color = "black") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = NULL,
         x = if (show_xlab) "Model" else NULL,
         y = y_label) +
    theme_minimal() +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

scatter_rgr_num <- make_scatter_rgr(
  df_rank_rgr_enriched %>% filter(data_type == "Numeric"), "Numeric")
scatter_rgr_cat <- make_scatter_rgr(
  df_rank_rgr_enriched %>% filter(data_type == "Categorical"), "Categorical")
scatter_rgr_all <- make_scatter_rgr(
  df_rank_rgr_enriched, "All", show_xlab = TRUE)

# ── 6c. Heatmap rank RGR by stratum ───────────────────────────────────────────

make_heatmap_rgr <- function(df_sub, show_xlab = FALSE) {
  df_sub <- df_sub %>%
    filter(model %in% spec_levels_rgr) %>%
    mutate(spec = as.character(model))
  
  df_perc <- tab_to_perc(df_sub, "spec", "RGR_rank") %>%
    mutate(spec = factor(spec, levels = rev(spec_levels_rgr)))
  
  ggplot(df_perc, aes(x = rank, y = spec, fill = perc)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "black",
                        name = "% Freq.", limits = c(0, 60)) +
    labs(title = NULL,
         x = if (show_xlab) "Rank" else NULL,
         y = NULL) +
    scale_y_discrete(limits = rev(spec_levels_rgr)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

heatmap_rgr_num <- make_heatmap_rgr(
  df_rank_rgr_enriched %>% filter(data_type == "Numeric"))
heatmap_rgr_cat <- make_heatmap_rgr(
  df_rank_rgr_enriched %>% filter(data_type == "Categorical"))
heatmap_rgr_all <- make_heatmap_rgr(
  df_rank_rgr_enriched, show_xlab = TRUE)

# ── 6d. RGR 2x3 panel ──────────────────────────────────────────────────────────

col_a_rgr <- (scatter_rgr_num / scatter_rgr_cat / scatter_rgr_all) +
  plot_layout(heights = c(1, 1, 1)) &
  theme(legend.position = "none")

col_b_rgr <- (heatmap_rgr_num / heatmap_rgr_cat / heatmap_rgr_all) +
  plot_layout(guides = "collect", heights = c(1, 1, 1)) &
  theme(legend.position = "right")

combined_rgr <- cowplot::plot_grid(
  col_a_rgr, col_b_rgr,
  labels     = c("a", "b"),
  label_size = 14,
  ncol       = 2,
  align      = "h",
  rel_widths = c(1, 1.08)
)

combined_rgr

# =============================================================================
# 7. FIGURE B — DELTA-RGR (improvement_RGR = mean_RGR - RGR_RF)
#    RF excluded (improvement_RGR(RF) = 0 by construction)
# =============================================================================

# ── 7a. Rank by delta-RGR (rank 1 = highest improvement) ─────────────────────

df_rank_imp <- df_rgr %>%
  filter(model %in% spec_levels_imp) %>%
  mutate(model = factor(model, levels = spec_levels_imp)) %>%
  group_by(dataset) %>%
  mutate(
    imp_RGR_rank = rank(-improvement_RGR, ties.method = "first")
  ) %>%
  ungroup()

df_rank_imp_enriched <- df_rank_imp %>%
  left_join(
    Info_248_df_safe %>%
      select(dataset, n_categorical_features, n_binary_features) %>%
      mutate(
        data_type_ = ifelse(
          n_categorical_features != 0 | n_binary_features != 0,
          "Categorical", "Numeric"
        )
      ),
    by = "dataset"
  )

# ── 7b. Scatter delta-RGR by stratum ──────────────────────────────────────────

make_scatter_imp <- function(df_sub, y_label, show_xlab = FALSE) {
  df_sub <- df_sub %>%
    filter(model %in% spec_levels_imp) %>%
    mutate(model = factor(model, levels = spec_levels_imp))
  
  ggplot(df_sub, aes(x = model, y = improvement_RGR)) +
    geom_jitter(
      data = df_sub %>%
        mutate(y_jitter = ifelse(
          improvement_RGR >  appiatimento_imp,  appiatimento_imp,
          ifelse(improvement_RGR < -appiatimento_imp, -appiatimento_imp,
                 improvement_RGR)
        )),
      aes(y = y_jitter),
      width = 0.2, height = 0,
      alpha = 0.5, size = 1.5, color = "darkgray"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    stat_summary(fun = mean, geom = "point",
                 shape = 4, size = 2, stroke = 1, color = "black") +
    coord_cartesian(ylim = c(-appiatimento_imp - 0.02, appiatimento_imp + 0.02)) +
    labs(title = NULL,
         x = if (show_xlab) "Model" else NULL,
         y = y_label) +
    theme_minimal() +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

scatter_imp_num <- make_scatter_imp(
  df_rank_imp_enriched %>% filter(data_type == "Numeric"), "Numeric")
scatter_imp_cat <- make_scatter_imp(
  df_rank_imp_enriched %>% filter(data_type == "Categorical"), "Categorical")
scatter_imp_all <- make_scatter_imp(
  df_rank_imp_enriched, "All", show_xlab = TRUE)

# ── 7c. Heatmap rank delta-RGR by stratum ─────────────────────────────────────

make_heatmap_imp <- function(df_sub, show_xlab = FALSE) {
  df_sub <- df_sub %>%
    filter(model %in% spec_levels_rgr) %>%
    mutate(spec = as.character(model))
  
  df_perc <- tab_to_perc(df_sub, "spec", "RGR_rank") %>%
    mutate(spec = factor(spec, levels = rev(spec_levels_rgr)))
  
  ggplot(df_perc, aes(x = rank, y = spec, fill = perc)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "black",
                        name = "% Freq.", limits = c(0, 60)) +
    labs(title = NULL,
         x = if (show_xlab) "Rank" else NULL,
         y = NULL) +
    scale_y_discrete(limits = rev(spec_levels_rgr)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

heatmap_imp_num <- make_heatmap_imp(
  df_rank_rgr_enriched %>% filter(data_type == "Numeric"))
heatmap_imp_cat <- make_heatmap_imp(
  df_rank_rgr_enriched %>% filter(data_type == "Categorical"))
heatmap_imp_all <- make_heatmap_imp(
  df_rank_rgr_enriched, show_xlab = TRUE)

# ── 7d. delta-RGR 2x3 panel ────────────────────────────────────────────────────

col_a_imp <- (scatter_imp_num  / scatter_imp_all) +
  plot_layout(heights = c(1, 1)) &
  theme(legend.position = "none")

col_b_imp <- (heatmap_imp_num  / heatmap_imp_all) +
  plot_layout(guides = "collect", heights = c(1, 1)) &
  theme(legend.position = "right")

combined_delta_rgr <- cowplot::plot_grid(
  col_a_imp, col_b_imp,
  labels     = c("a", "b"),
  label_size = 14,
  ncol       = 2,
  align      = "h",
  rel_widths = c(1, 1.08)
)

combined_delta_rgr


ggplot2::ggsave(
  filename = file.path(FIGDIR, "fig_RGR_delta_2x2.pdf"),
  plot     = combined_delta_rgr,
  device   = cairo_pdf,
  width    = 170, height = 133, units = "mm", dpi = 200
)

# =============================================================================
# 8. SAVE OUTPUT
# =============================================================================

# dir.create("output_analysis", showWarnings = FALSE)

# ggplot2::ggsave(
#   filename = "output_analysis/fig_RGR_2x3.pdf",
#   plot     = combined_rgr,
#   device   = cairo_pdf,
#   width    = 170, height = 200, units = "mm", dpi = 200
# )
#
# ggplot2::ggsave(
#   filename = "output_analysis/fig_deltaRGR_2x3.pdf",
#   plot     = combined_delta_rgr,
#   device   = cairo_pdf,
#   width    = 170, height = 200, units = "mm", dpi = 200
# )

saveRDS(tabella_RGR, file.path(TABDIR, "tabella_RGR_imp_benchmark.rds"))
cat("\nEND fig5_benchmark_rgr.R\n")
