# ==============================================================================
# Shared helpers for the run_benchmark_<MODEL>.R "from scratch" generation
# scripts (competing models trained on the 248 real datasets).
# Sourced by each run_benchmark_*.R script; not meant to be run directly.
# ==============================================================================

sanitize_colnames <- function(nomi) {
  nomi <- gsub("[><=]", "_", nomi)
  nomi <- gsub("[^[:alnum:]_]", "_", nomi)
  nomi <- gsub("^[0-9]", "X", nomi)
  nomi <- gsub("weight", "WEIGHT", nomi)
  nomi <- gsub("time",   "TIME",   nomi)
  nomi <- make.names(nomi, unique = TRUE)
  return(nomi)
}

rename_except_target <- function(df, target_name = "target") {
  if (!(target_name %in% colnames(df))) {
    warning("Target column not found. Proceeding without exclusion.")
    target_name <- NULL
  }
  cols_to_rename <- setdiff(colnames(df), target_name)
  generate_names <- function(n) {
    base <- LETTERS
    if (n <= length(base)) return(paste0("var", base[1:n]))
    comb <- c(base, outer(base, base, paste0))
    return(paste0("var", comb[1:n]))
  }
  new_names <- generate_names(length(cols_to_rename))
  colnames(df)[colnames(df) %in% cols_to_rename] <- new_names
  return(df)
}

# Rank Graduation Accuracy (RGA), rank-based measure in [0, 1]:
#   1   = perfect concordance (predicted ranking = true ranking)
#   0.5 = ~random behaviour
#   0   = perfect discordance (predicted ranking = reverse of the true one)
RGA_metric <- function(y_true, y_pred) {
  y_true <- as.numeric(y_true)
  y_pred <- as.numeric(y_pred)
  ok     <- is.finite(y_true) & is.finite(y_pred)
  y_true <- y_true[ok]
  y_pred <- y_pred[ok]
  n <- length(y_true)
  if (n < 2) return(NA_real_)
  r_true     <- order(y_true)
  r_true_rev <- order(y_true, decreasing = TRUE)
  r_pred     <- order(y_pred)
  S_pred  <- sum((1:n) * y_true[r_pred])
  S_best  <- sum((1:n) * y_true[r_true])
  S_worst <- sum((1:n) * y_true[r_true_rev])
  num <- S_pred - S_worst
  den <- S_best - S_worst
  if (den == 0) return(NA_real_)
  num / den
}

# Deterministic perturbation (5%-95% swap), as in Babaei et al. 2025:
# retrain on this to get the "perturbed" robustness cycle. Variables with
# <= 10 unique values (categorical/low-cardinality) are left untouched.
perturb_variables <- function(x_to_perturb, perturbation_percentage = 0.05) {
  x_pert  <- as.data.frame(x_to_perturb)
  n_total <- nrow(x_pert)
  for (j in seq_len(ncol(x_pert))) {
    col <- as.numeric(x_pert[, j])
    if (length(unique(col)) <= 10) next
    sorted_idx <- order(col)
    p5_n    <- as.integer(ceiling(perturbation_percentage * n_total))
    p95_n   <- as.integer(ceiling((1 - perturbation_percentage) * n_total))
    n_upper <- n_total - p95_n
    n       <- min(p5_n, n_upper)
    if (n == 0) next
    lower_idx <- sorted_idx[seq_len(n)]
    upper_idx <- rev(sorted_idx[(n_total - n + 1):n_total])
    col_new            <- col
    col_new[lower_idx] <- col[upper_idx]
    col_new[upper_idx] <- col[lower_idx]
    x_pert[, j]        <- col_new
  }
  x_pert
}

# Standard preprocessing pipeline shared by every model:
# clean column names -> dummy-encode factors/characters -> sanitize names
# -> (optionally) rename predictors to generic varA/varB/... .
preprocess_benchmark_df <- function(df, do_rename = TRUE) {
  colnames(df) <- gsub("-", "_", colnames(df))
  colnames(df) <- gsub(" ", "_", colnames(df))
  if (any(sapply(df, is.factor) | sapply(df, is.character))) {
    df <- fastDummies::dummy_cols(df, remove_selected_columns = TRUE,
                                  remove_first_dummy = TRUE)
  }
  colnames(df) <- sanitize_colnames(colnames(df))
  if (do_rename) df <- rename_except_target(df)
  df
}

# 70/30 train/test partitions: 200 pre-generated splits via
# caret::createDataPartition (set.seed(10)), or 3500-row random subsamples
# for datasets with more than 5000 rows. Repetitions 1:50 index into this
# matrix (column = repetition).
make_partition_matrix <- function(df) {
  if (nrow(df) > 5000) {
    partition_matrix <- matrix(NA, nrow = 3500, ncol = 200)
    set.seed(10)
    for (ii in 1:200) partition_matrix[, ii] <- sample(nrow(df), 3500)
  } else {
    set.seed(10)
    partition_matrix <- caret::createDataPartition(df$target, p = 0.7,
                                                   times = 200, list = FALSE)
  }
  partition_matrix
}
