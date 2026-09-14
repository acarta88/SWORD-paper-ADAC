# =============================================================================
# Fig. 6 (fig:rge_vs_sword) — RGE figure for the paper -> fig_RGE_vs_SWORD_v3.pdf
# Figure stage only: reads output_RGE/df_RGE_wide.rds (aggregated, ~0.3 MB).
# No recomputation from the raw RGE predictions (not deposited in this
# repository -- see run_rge_predictions.R).
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
PROJ    <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
DATA    <- file.path(PROJ, "data")
FIGDIR  <- file.path(PROJ, "figures")
OUT_RGE <- file.path(SCRIPT_DIR, "output_RGE")
if (!dir.exists(FIGDIR)) dir.create(FIGDIR, recursive = TRUE)
if (!interactive()) pdf(NULL)

library(dplyr)
library(ggplot2)
library(tidyr)

out_dir <- OUT_RGE
load(file.path(DATA, "benchmark_248_dataset_list_ADAC.rda"))
# =============================================================================
# 1. Load data
# =============================================================================

df_wide <- readRDS(file.path(out_dir, "df_RGE_wide.rds"))

#eliminate df with only 2 predictors
ds_3var <- names(Filter(function(df) ncol(df) == 3, benchmark_248_dataset_list))
ds_3var
df_wide <- df_wide[!(df_wide$dataset %in% ds_3var), ]

rge_cols <- grep("^RGE_", names(df_wide), value = TRUE)
df_rge   <- df_wide[, rge_cols]
names(df_rge) <- sub("^RGE_", "", names(df_rge))

cat("Available models:", paste(names(df_rge), collapse = ", "), "\n")

# =============================================================================
# 2. Generic function: 5x2 facet of "others" vs "reference"
# =============================================================================

make_paper_figure <- function(df, ref_model, out_file,
                              x_lab = NULL, y_lab = NULL,
                              ncols = 5) {
  
  others <- setdiff(names(df), ref_model)
  
  df_long <- df %>%
    rename(REF = all_of(ref_model)) %>%
    pivot_longer(cols = all_of(others),
                 names_to  = "model",
                 values_to = "RGE_other") %>%
    filter(is.finite(REF), is.finite(RGE_other))
  
  labels <- df_long %>%
    group_by(model) %>%
    summarise(
      beta = coef(lm(RGE_other ~ REF))[2],
      r2   = summary(lm(RGE_other ~ REF))$r.squared,
      n    = n(),
      .groups = "drop"
    ) %>%
    mutate(label = sprintf("atop(beta == %.2f, R^2 == %.2f)", beta, r2))
  
  panel_order <- c("SWORD", "ODRF", "EXTree", "RF", "RRF",
                   "aorsfNET", "aorsf", "SPORF", "CF", "ROT", "SVMlin")
  panel_order <- panel_order[panel_order %in% others]
  df_long$model <- factor(df_long$model, levels = panel_order)
  labels$model  <- factor(labels$model,  levels = panel_order)
  
  x_label <- if (!is.null(x_lab)) x_lab else paste("RGE —", ref_model)
  y_label <- if (!is.null(y_lab)) y_lab else "RGE — other model"
  
  p <- ggplot(df_long, aes(x = REF, y = RGE_other)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "grey70", linewidth = 0.4) +
    geom_point(alpha = 0.20, size = 0.7, color = "black") +
    geom_smooth(method = "lm", se = FALSE,
                color     = "grey30",
                #linetype  = "dashed",
                linewidth = 0.9,
                fullrange = TRUE) +
    geom_text(data = labels,
              aes(x = 0.03, y = 0.96, label = label),
              hjust = 0, vjust = 1,
              size  = 3.5,
              color = "black",
              parse = TRUE,
              inherit.aes = FALSE) +
    facet_wrap(~ model, ncol = ncols) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       labels = c("0", "0.5", "1"),
                       expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       labels = c("0", "0.5", "1"),
                       expand = c(0.01, 0.01)) +
    labs(x = x_label,
         y = y_label) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_blank(),
      strip.text       = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "grey92", color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      axis.title       = element_text(size = 10),
      axis.text        = element_text(size = 8),
      plot.margin      = margin(4, 4, 4, 4)
    )
  
  ggsave(
    filename = file.path(FIGDIR, out_file),
    plot     = p,
    width    = 10,
    height   = 5,
    units    = "in",
    device   = "pdf"
  )
  
  cat(sprintf("Saved: %s  (10 x 5 in, %d panels)\n", out_file, length(others)))
}

# =============================================================================
# 3. Fig 1 – all vs RF
# =============================================================================

make_paper_figure(
  df        = df_rge,
  ref_model = "RF",
  out_file  = "fig_RGE_vs_RF_v2.pdf",
  x_lab     = "RGE  —  Random Forests",
  y_lab     = "RGE  —  Other Methods"
)

# =============================================================================
# 4. Fig 2 – all vs SWORD
# =============================================================================

make_paper_figure(
  df        = df_rge,
  ref_model = "SWORD",
  out_file  = "fig_RGE_vs_SWORD_v3.pdf",
  x_lab     = "RGE  —  SWORD",
  y_lab     = "RGE  —  Other Methods"
)

cat("\nDone.\n")