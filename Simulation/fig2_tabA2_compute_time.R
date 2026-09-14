# =============================================================================
# Fig. 2 (fig:compute_time) + Table A2 (tab:compute_time_scaling)
# Computational time per tree — W_scheme = scale
# Config: Pearson | Random | C=1 | rf_var=100%
# n_obs x features grid (3x3)
#
# Output: figures/fig_compute_time_paper_72.pdf
#         (heatmap without SD + log-log, side by side; same original style)
# =============================================================================

library(dplyr)
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
CACHE   <- file.path(PROJ, "results_cache")
dir_out <- file.path(PROJ, "figures")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)
if (!interactive()) pdf(NULL)

# ── Load only scale ────────────────────────────────────────────────────────

path_cache <- file.path(CACHE, "df_SWORD72_cached.rds")



ricostruisci_df_old <- function(cartella) {
  file_rds <- list.files(cartella, pattern = "\\.rds$", full.names = TRUE)
  file_rds <- file_rds[basename(file_rds) != "risultati_completi.rds"]
  cat("  File .rds in", basename(cartella), ":", length(file_rds), "\n")
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
  cat("  File .rds in", basename(cartella), ":", length(file_rds), "\n")
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
df_finale_SWORD_scale=df

df <- df_finale_SWORD_scale %>%
  mutate(
    features    = as.integer(features),
    n_obs       = as.integer(n_obs),
    rf_var      = as.integer(rf_var),
    relation    = factor(relation,    levels = c("Pearson", "MI")),
    Cost_C      = factor(as.character(Cost_C), levels = c("0.1", "1", "10")),
    rf_var_frac = factor(
      cut(rf_var / features, breaks = c(0, 0.45, 0.75, Inf),
          labels = c("30%", "60%", "100%"), right = TRUE),
      levels = c("30%", "60%", "100%")
    ),
    rand_label  = factor(ifelse(rand_ntopcor, "Random", "Fixed"),
                         levels = c("Fixed", "Random"))
  )

# Target config
df_target <- df %>%
  filter(
    W_scheme == "scale",
    relation    == "Pearson",
    rand_label  == "Random",
    Cost_C      == "1",
    rf_var_frac == "100%",
    !is.na(time_tree_mean)
  )

cat(sprintf("Rows (Pe | Random | C=1 | 100%%): %d\n", nrow(df_target)))

# Aggregation: mean over noise, nonlin, cat, seed -> n_obs x features grid
tempo_grid <- df_target %>%
  group_by(n_obs, features) %>%
  summarise(
    tempo_medio = mean(time_tree_mean, na.rm = TRUE),
    tempo_sd    = sd(time_tree_mean,   na.rm = TRUE),
    n_dgp       = n(),
    .groups = "drop"
  )

cat("\n=== Mean tree-building time (seconds) — n_obs x features grid ===\n")
print(tempo_grid %>% arrange(n_obs, features), n = 9)

# Scaling ratios (n=1000 vs n=100)
cat("\n=== Scaling ratios: n=1000 vs n=100 ===\n")
tempo_grid %>%
  group_by(features) %>%
  summarise(
    t_100  = tempo_medio[n_obs == 100],
    t_1000 = tempo_medio[n_obs == 1000],
    ratio  = t_1000 / t_100,
    .groups = "drop"
  ) %>%
  print()

# ── FIGURE A: Heatmap n_obs x features (without SD) ────────────────────────────

mid_g <- median(tempo_grid$tempo_medio)

p_heat <- ggplot(
  tempo_grid,
  aes(x = factor(n_obs), y = factor(features), fill = tempo_medio)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", tempo_medio)), size = 4) +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "#F7F7F7",
    high     = "#D6604D",
    midpoint = mid_g,
    name     = "Time (s)"
  ) +
  labs(
    x = "Sample size",
    y = "Number of predictors"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid      = element_blank(),
    legend.position = "right"
  )

# ── FIGURE B: Log-log plot with reference lines ──────────────────────────

# Lines anchored at the point (n=100, p=5)
ref_base <- tempo_grid %>%
  filter(features == 5, n_obs == 100) %>%
  pull(tempo_medio)

ref_lines <- data.frame(
  n_obs = rep(c(100, 500, 1000), 2),
  tipo  = rep(c("Linear (\u03b1 = 1)", "Quadratic (\u03b1 = 2)"), each = 3)
) %>%
  mutate(
    tempo = ifelse(
      tipo == "Linear (\u03b1 = 1)",
      ref_base * (n_obs / 100)^1,
      ref_base * (n_obs / 100)^2
    )
  )

p_loglog <- ggplot(
  tempo_grid,
  aes(
    x     = n_obs,
    y     = tempo_medio,
    color = factor(features),
    shape = factor(features),
    group = factor(features)
  )
) +
  # ── Reference lines (no legend) ────────────────────────────
  geom_line(
    data        = ref_lines %>% filter(tipo == "Linear (\u03b1 = 1)"),
    aes(x = n_obs, y = tempo),
    linetype    = "dashed",
    color       = "grey55",
    linewidth   = 0.65,
    inherit.aes = FALSE
  ) +
  geom_line(
    data        = ref_lines %>% filter(tipo == "Quadratic (\u03b1 = 2)"),
    aes(x = n_obs, y = tempo),
    linetype    = "dotted",
    color       = "grey55",
    linewidth   = 0.65,
    inherit.aes = FALSE
  ) +

  # ── Actual curves ──────────────────────────────────────────────────────
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +

  # ── Scale log-log ─────────────────────────────────────────────────────
  scale_x_log10(
    breaks = c(100, 500, 1000),
    labels = c("100", "500", "1000")
  ) +
  scale_y_log10() +

  # ── Colors + shape ────────────────────────────────────────────────────
  scale_color_manual(
    values = c("#1a6faf", "#e07020", "#228b22"),
    name   = "Number of predictors"
  ) +
  scale_shape_manual(
    values = c(16, 17, 15),
    name   = "Number of predictors"
  ) +

  annotation_logticks(sides = "bl", size = 0.3) +

  labs(
    x = "Sample size",
    y = "Mean tree time"
  ) +

  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor     = element_blank(),
    # legend inside
    legend.position      = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = "white", color = "grey70"),
    legend.key           = element_rect(fill = "white", color = NA),
    legend.title         = element_text(size = 11),
    legend.text          = element_text(size = 10)
  ) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1.1))
  )

# ── Combine and save ───────────────────────────────────────────────────────────

fig_combined <- ggarrange(
  p_heat, p_loglog,
  ncol       = 2,
  widths     = c(1, 1.2),
  labels     = c("a", "b"),
  font.label = list(size = 14, face = "bold"),
  label.x    = 0.02,
  label.y    = 0.98
)

out_path <- file.path(dir_out, "fig_compute_time_paper_72.pdf")
ggsave(out_path, fig_combined, width = 10, height = 4)

cat(sprintf("\nSaved: %s\n", out_path))
fig_combined

# =============================================================================
# LaTeX appendix table: scaling exponents kappa (in n) and beta (in m)
# =============================================================================

m_vals <- sort(unique(tempo_grid$features))  # c(5, 10, 50)
n_vals <- sort(unique(tempo_grid$n_obs))      # c(100, 500, 1000)

# --- kappa: scaling in n, for each level of m ---
kappa_tab <- tempo_grid %>%
  arrange(features, n_obs) %>%
  group_by(features) %>%
  summarise(
    t_100  = tempo_medio[n_obs == 100],
    t_500  = tempo_medio[n_obs == 500],
    t_1000 = tempo_medio[n_obs == 1000],
    ratio_n = t_1000 / t_100,
    kappa   = coef(lm(log(c(t_100, t_500, t_1000)) ~
                      log(c(100, 500, 1000))))[2],
    .groups = "drop"
  )

# --- beta: scaling in m, for each level of n ---
beta_tab <- tempo_grid %>%
  arrange(n_obs, features) %>%
  group_by(n_obs) %>%
  summarise(
    t_5  = tempo_medio[features == 5],
    t_10 = tempo_medio[features == 10],
    t_50 = tempo_medio[features == 50],
    ratio_m = t_50 / t_5,
    beta    = coef(lm(log(c(t_5, t_10, t_50)) ~
                      log(c(5, 10, 50))))[2],
    .groups = "drop"
  )

cat("\n=== Scaling in n (kappa) by level of m ===\n")
print(kappa_tab %>% mutate(across(where(is.numeric), ~ round(., 3))))

cat("\n=== Scaling in m (beta) by level of n ===\n")
print(beta_tab %>% mutate(across(where(is.numeric), ~ round(., 3))))

