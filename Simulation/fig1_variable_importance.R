# ==============================================================================
# Fig. 1 (fig:Var_imp) — VIPLOT: Scatter VI SWORD vs RF — seed_1
# Recommended SWORD config: W_std (scale) | Pearson | Random gamma | C=1 | alpha=100%
# Seed: only seed_1
# Normalization: VI / sum(VI) per dataset -> values in [0, 1], sum = 1
# Output: figures/VIPLOT_sword_rf_seed1_puliti_scale.pdf
# ==============================================================================

library(dplyr)
library(ggplot2)

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
PROJ      <- normalizePath(file.path(SCRIPT_DIR, ".."))
SCALE_DIR <- file.path(SCRIPT_DIR, "risultati_full_simulation_scale")
RF_DIR    <- file.path(SCRIPT_DIR, "risultati_RF_simulation")
dir_out   <- file.path(PROJ, "figures")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)
if (!interactive()) pdf(NULL)


# ==============================================================================
# 1. LOAD AND FILTER VI SWORD
#    W_scheme=scale | Pearson | rand_ntopcor=TRUE | Cost_C=1 | rf_var_frac=100%
#    New structure: obj$combos -> list of elements with x$res and x$VI
#    File: scale__seed_1__<n_obs>_n_obs__<feat>_n_feat__...
#    val: [NA, seed, n_obs, features, noise_prop, nonlin, cat_prop, error_sd]
# ==============================================================================

sword_vi_cache <- file.path(SCALE_DIR, "VI_SWORD_seed1_scale_cached.rds")

if (file.exists(sword_vi_cache)) {

  cat("Loading SWORD VI (seed_1) from cache...\n")
  vi_sword_raw <- readRDS(sword_vi_cache)

} else {

  cat("SWORD VI cache not found - rebuilding from .rds (needs the per-dataset raw files)\n")

  file_rds_sword <- list.files(SCALE_DIR, pattern = "\\.rds$",
                                full.names = TRUE)
  file_rds_sword <- file_rds_sword[basename(file_rds_sword) != "risultati_completi.rds"]
  file_rds_sword <- file_rds_sword[grepl("^scale__seed_1__", basename(file_rds_sword))]
  cat("  Files found:", length(file_rds_sword), "\n")

  vi_sword_raw <- do.call(rbind, Filter(Negate(is.null), lapply(file_rds_sword, function(f) {
    nm    <- gsub("\\.rds$", "", basename(f))
    el    <- strsplit(nm, "__")[[1]]
    val   <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", el)))
    # val: [NA(scale), seed, n_obs, features, noise_prop, nonlin, cat_prop, error_sd]
    n_feat <- val[4]

    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(obj)) return(NULL)

    # New structure: obj$combos is the list of combos
    combo_list <- if (!is.null(obj$combos)) obj$combos else
      Filter(function(x) is.list(x) && !is.null(x$res), obj)
    if (length(combo_list) == 0) return(NULL)

    do.call(rbind, Filter(Negate(is.null), lapply(combo_list, function(x) {
      rel   <- x$res$relation
      rnd   <- as.logical(x$res$rand_ntopcor)
      costC <- as.numeric(as.character(x$res$Cost_C))
      frac  <- x$res$rf_var / n_feat

      # Filter: Pearson | Random | C=1 | 100%
      if (!isTRUE(rel == "Pearson")) return(NULL)
      if (!isTRUE(rnd == TRUE))     return(NULL)
      if (!isTRUE(costC == 1))      return(NULL)
      if (frac < 0.75)               return(NULL)

      vi <- x$VI
      if (is.null(vi) || length(vi) == 0 || is.null(names(vi))) return(NULL)

      data.frame(
        variable   = names(vi),
        VI_SWORD   = as.numeric(vi),
        n_obs      = val[3],
        features   = n_feat,
        noise_prop = val[5],
        nonlin     = val[6],
        cat_prop   = val[7],
        error_sd   = val[8],
        stringsAsFactors = FALSE
      )
    })))
  })))

}

cat("  VI SWORD rows extracted:", nrow(vi_sword_raw), "\n")
cat("  Unique datasets SWORD:", length(unique(with(vi_sword_raw,
    paste(n_obs, features, noise_prop, nonlin, cat_prop, error_sd)))), "\n")

# Normalize VI SWORD per dataset (divide by the sum -> sums to 1)
vi_sword_raw <- vi_sword_raw %>%
  group_by(n_obs, features, noise_prop, nonlin, cat_prop, error_sd) %>%
  mutate(VI_SWORD = {
    s <- sum(VI_SWORD, na.rm = TRUE)
    if (s > 1e-12) VI_SWORD / s else VI_SWORD
  }) %>%
  ungroup()


# ==============================================================================
# 2. LOAD VI RF (seed_1)
# ==============================================================================

rf_vi_cache <- file.path(RF_DIR, "VI_RF_seed1_cached.rds")

if (file.exists(rf_vi_cache)) {

  cat("\nLoading RF VI (seed_1) from cache...\n")
  vi_rf_raw <- readRDS(rf_vi_cache)

} else {

  cat("\nRF VI cache not found - rebuilding from .rds (needs the per-dataset raw files)\n")

  file_rds_rf <- list.files(RF_DIR, pattern = "\\.rds$",
                             full.names = TRUE)
  file_rds_rf <- file_rds_rf[basename(file_rds_rf) != "risultati_completi.rds"]
  file_rds_rf <- file_rds_rf[grepl("^seed_1__", basename(file_rds_rf))]
  cat("  Files found:", length(file_rds_rf), "\n")

  vi_rf_raw <- do.call(rbind, Filter(Negate(is.null), lapply(file_rds_rf, function(f) {
    nm  <- gsub("\\.rds$", "", basename(f))
    el  <- strsplit(nm, "__")[[1]]
    val <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", el)))

    res <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(res)) return(NULL)

    # Take the first non-null element (RF has a single config per file)
    do.call(rbind, Filter(Negate(is.null), lapply(res, function(x) {
      vi <- x$VI
      if (is.null(vi)) return(NULL)
      if (is.matrix(vi) || is.data.frame(vi)) {
        vi_names <- rownames(vi)
        vi_vals  <- as.numeric(vi[, 1])
      } else {
        vi_names <- names(vi)
        vi_vals  <- as.numeric(vi)
      }
      if (is.null(vi_names) || length(vi_names) == 0) return(NULL)
      # RF MDI may have negative values -> clip to 0
      vi_vals <- pmax(vi_vals, 0)

      data.frame(
        variable   = vi_names,
        VI_RF      = vi_vals,
        n_obs      = val[2],
        features   = val[3],
        noise_prop = val[4],
        nonlin     = val[5],
        cat_prop   = val[6],
        error_sd   = val[7],
        stringsAsFactors = FALSE
      )
    })))
  })))

}

cat("  VI RF rows extracted:", nrow(vi_rf_raw), "\n")
cat("  Unique datasets RF:", length(unique(with(vi_rf_raw,
    paste(n_obs, features, noise_prop, nonlin, cat_prop, error_sd)))), "\n")

# Normalize VI RF per dataset
vi_rf_raw <- vi_rf_raw %>%
  group_by(n_obs, features, noise_prop, nonlin, cat_prop, error_sd) %>%
  mutate(VI_RF = {
    s <- sum(VI_RF, na.rm = TRUE)
    if (s > 1e-12) VI_RF / s else VI_RF
  }) %>%
  ungroup()


# ==============================================================================
# 3. MERGE: build df_VI_combined
# ==============================================================================

id_cols <- c("variable", "n_obs", "features", "noise_prop",
             "nonlin", "cat_prop", "error_sd")

df_VI_combined <- merge(
  vi_rf_raw[,   c(id_cols, "VI_RF")],
  vi_sword_raw[, c(id_cols, "VI_SWORD")],
  by = id_cols
)

cat("\ndf_VI_combined: ", nrow(df_VI_combined), "rows,",
    length(unique(with(df_VI_combined,
      paste(n_obs, features, noise_prop, nonlin, cat_prop, error_sd)))),
    "dataset\n")

# Tag variables
df_VI_combined <- df_VI_combined %>%
  mutate(
    noise_var = grepl("noise", variable),
    cat_var   = grepl("[Cc]at",  variable),
    NL_var    = grepl("funz|nl|NL", variable)
  )

cat("Noise variables:", sum(df_VI_combined$noise_var),
    "| informative:", sum(!df_VI_combined$noise_var), "\n")


# ==============================================================================
# 4. REGRESSION MODELS
# ==============================================================================

x_start <- 0
x_end   <- 1

regr_line <- function(model, x_vals) {
  predict(model, newdata = data.frame(VI_RF = x_vals))
}

mod_total <- lm(VI_SWORD ~ VI_RF, data = df_VI_combined)
mod_noise <- lm(VI_SWORD ~ VI_RF, data = subset(df_VI_combined,  noise_var))
mod_info  <- lm(VI_SWORD ~ VI_RF, data = subset(df_VI_combined, !noise_var))

cat("\nTotal model    — R²=", round(summary(mod_total)$r.squared, 4),
    " slope=", round(coef(mod_total)[2], 4), "\n")
cat("Noise model     — R²=", round(summary(mod_noise)$r.squared, 4),
    " slope=", round(coef(mod_noise)[2], 4), "\n")
cat("Info model      — R²=", round(summary(mod_info)$r.squared,  4),
    " slope=", round(coef(mod_info)[2],  4), "\n")


# ==============================================================================
# 5. FIGURE: scatter + bisector + regression lines
#    Reproduces exactly the structure of VIPLOT_sword_rf_seed1_puliti
# ==============================================================================

VI_plot <- ggplot(df_VI_combined,
                  aes(x = VI_RF, y = VI_SWORD, color = noise_var)) +

  # Points
  geom_point(alpha = 0.5, size = 1.5) +

  # Bisector y = x (dotted black)
  annotate("segment",
    x = x_start, xend = x_end, y = x_start, yend = x_end,
    linetype = "dotted", color = "black", linewidth = 0.7
  ) +

  # Global regression line (solid black)
  annotate("segment",
    x = x_start, xend = x_end,
    y = regr_line(mod_total, x_start), yend = regr_line(mod_total, x_end),
    color = "black", linewidth = 1
  ) +

  # Noise regression line (dashed red)
  annotate("segment",
    x = x_start, xend = x_end,
    y = regr_line(mod_noise, x_start), yend = regr_line(mod_noise, x_end),
    color = "red", linetype = "dashed", linewidth = 0.9
  ) +

  # Informative regression line (dashed blue)
  annotate("segment",
    x = x_start, xend = x_end,
    y = regr_line(mod_info, x_start), yend = regr_line(mod_info, x_end),
    color = "blue", linetype = "dashed", linewidth = 0.9
  ) +

  # Colors: red = noise, blue = informative
  scale_color_manual(
    values = c(
      "TRUE"  = rgb(1, 0, 0, 0.6),
      "FALSE" = rgb(0, 0, 1, 0.6)
    ),
    labels = c(
      "FALSE" = "Informative variable",
      "TRUE"  = "Noise variable"
    ),
    name = "Variable type"
  ) +

  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  coord_fixed(ratio = 1) +

  labs(
    x = "RF Variable Importance",
    y = "SWORD Variable Importance"
  ) +

  theme_minimal(base_size = 16) +
  theme(
    legend.position = "top",
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 13)
  )

# Save
ggsave(
  file.path(dir_out, "VIPLOT_sword_rf_seed1_puliti_scale.pdf"),
  VI_plot,
  width = 6, height = 6, units = "in", dpi = 600
)
cat("\nSaved: figures/VIPLOT_sword_rf_seed1_puliti_scale.pdf\n")

# Also show in RStudio
VI_plot


# ==============================================================================
# 6. GINI COEFFICIENT: VI concentration - OIW-VI vs RF MDI
# ==============================================================================

gini_coef <- function(x) {
  x <- sort(x[!is.na(x) & x >= 0])
  n <- length(x); s <- sum(x)
  if (n == 0 || s == 0) return(NA_real_)
  2 * sum(seq_len(n) * x) / (n * s) - (n + 1) / n
}

gini_per_ds <- df_VI_combined %>%
  group_by(n_obs, features, noise_prop, nonlin, cat_prop, error_sd) %>%
  summarise(
    gini_sword = gini_coef(VI_SWORD),
    gini_rf    = gini_coef(VI_RF),
    .groups = "drop"
  ) %>%
  filter(!is.na(gini_sword), !is.na(gini_rf))

cat("\n=== Gini coefficient: OIW-VI vs RF MDI ===\n")
cat(sprintf("N dataset: %d\n",                       nrow(gini_per_ds)))
cat(sprintf("Average Gini OIW-VI : %.3f\n",            mean(gini_per_ds$gini_sword)))
cat(sprintf("Average Gini RF MDI : %.3f\n",            mean(gini_per_ds$gini_rf)))
cat(sprintf("Mean diff (OIW-VI - RF MDI): %.3f\n",  mean(gini_per_ds$gini_sword - gini_per_ds$gini_rf)))

wt_gini <- wilcox.test(gini_per_ds$gini_sword, gini_per_ds$gini_rf, paired = TRUE)
cat(sprintf("Wilcoxon signed-rank: W=%.0f, p=%.4f\n",
            wt_gini$statistic, wt_gini$p.value))


# ==============================================================================
# 7. GROUND TRUTH VALIDATION - mean_VI_info vs mean_VI_noise for OIW-VI
#    Only datasets with variables of both categories (noise_prop > 0)
#    H1: mean_VI_info > mean_VI_noise within OIW-VI
#    -> OIW-VI correctly identifies the informative variables
# ==============================================================================

sep <- df_VI_combined %>%
  filter(noise_prop > 0) %>%
  group_by(n_obs, features, noise_prop, nonlin, cat_prop, error_sd) %>%
  summarise(
    mean_vi_info  = mean(VI_SWORD[!noise_var], na.rm = TRUE),
    mean_vi_noise = mean(VI_SWORD[ noise_var], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(mean_vi_info), !is.na(mean_vi_noise))

cat("\n=== Ground truth validation: mean_VI_info vs mean_VI_noise (OIW-VI) ===\n")
cat(sprintf("Datasets analyzed (noise_prop > 0): %d\n", nrow(sep)))
cat(sprintf("%% datasets with mean_VI_info > mean_VI_noise: %.1f%%\n",
            100 * mean(sep$mean_vi_info > sep$mean_vi_noise)))
cat(sprintf("Median mean_VI_info  : %.5f\n", median(sep$mean_vi_info)))
cat(sprintf("Median mean_vi_noise : %.5f\n", median(sep$mean_vi_noise)))

wt_sep <- wilcox.test(sep$mean_vi_info, sep$mean_vi_noise,
                       paired = TRUE, alternative = "greater")
cat(sprintf("Wilcoxon signed-rank (info > noise): V=%.0f, p=%.2e\n",
            wt_sep$statistic, wt_sep$p.value))

# Effect size r_rb
diffs_sep <- sep$mean_vi_info - sep$mean_vi_noise
r_sep     <- rank(abs(diffs_sep))
W_plus    <- sum(r_sep[diffs_sep > 0])
W_minus   <- sum(r_sep[diffs_sep < 0])
r_rb_sep  <- (W_plus - W_minus) / (W_plus + W_minus)
cat(sprintf("Effect size r_rb = %.3f\n", r_rb_sep))







