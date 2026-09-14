# =============================================================================
# Table 4 (tab:summary_realds) — Overview table of real datasets
# Cross-tab: N. features (4 quartile-based classes) x N. instances (4 quartile-based classes)
# for the 248 benchmark datasets. Input: data/Info_df_ADAC.rda (Info_248_df).
# Output: tables/table_summary_realds.csv  +  tables/table_summary_realds.tex
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
PROJ   <- normalizePath(file.path(SCRIPT_DIR, ".."))
DATA   <- file.path(PROJ, "data")
TABDIR <- file.path(PROJ, "tables")
if (!dir.exists(TABDIR)) dir.create(TABDIR, recursive = TRUE)

load(file.path(DATA, "Info_df_ADAC.rda"))   # -> Info_248_df

d <- Info_248_df
stopifnot(all(c("n_features", "n_instances") %in% names(d)))

# --- Quartile-based classes (internal: 25/50/75%) ---
make_classes <- function(x) {
  qs <- quantile(x, probs = c(0.25, 0.50, 0.75), type = 7)
  brks <- c(min(x), floor(qs), max(x))
  brks <- unique(brks)
  # "lo-hi" labels on groups
  cut(x, breaks = c(-Inf, floor(qs), Inf), right = TRUE,
      labels = c(
        sprintf("%d–%d", min(x),            floor(qs)[1]),
        sprintf("%d–%d", floor(qs)[1] + 1,  floor(qs)[2]),
        sprintf("%d–%d", floor(qs)[2] + 1,  floor(qs)[3]),
        sprintf("%d–%d", floor(qs)[3] + 1,  max(x))
      ))
}

d$feat_class <- make_classes(d$n_features)
d$inst_class <- make_classes(d$n_instances)

# --- Cross-tab with totals ---
tab <- table(Features = d$feat_class, Instances = d$inst_class)
tab_tot <- addmargins(tab)

cat("\n=== tab:summary_realds — cross-tab (rows=Features, columns=Instances) ===\n")
print(tab_tot)
cat(sprintf("\nTotal datasets: %d\n", nrow(d)))

# --- Save outputs ---
df_out <- as.data.frame.matrix(tab_tot)
write.csv(df_out, file.path(TABDIR, "table_summary_realds.csv"))

# Minimal LaTeX (table body)
sink(file.path(TABDIR, "table_summary_realds.tex"))
cat("% tab:summary_realds — script: tab4_summary_realds.R\n")
inst_levels <- colnames(tab)
cat(paste0("Features \\textbackslash Instances & ",
           paste(inst_levels, collapse = " & "), " & Total \\\\\n"))
for (fl in rownames(tab)) {
  row_counts <- tab[fl, ]
  cat(paste0(fl, " & ", paste(row_counts, collapse = " & "),
             " & ", sum(row_counts), " \\\\\n"))
}
cat(paste0("Total & ", paste(colSums(tab), collapse = " & "),
           " & ", sum(tab), " \\\\\n"))
sink()

cat(sprintf("\nSaved in %s:\n  - table_summary_realds.csv\n  - table_summary_realds.tex\n", TABDIR))
