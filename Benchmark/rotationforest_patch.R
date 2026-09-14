# ==============================================================================
# Rotation Forest (Rodriguez et al. 2006) -- from-scratch implementation.
#
# Rotation Forest has no maintained CRAN package for this regression use
# case, so this is a from-scratch implementation built on rpart. Sourced by
# both run_benchmark_ROT.R and run_rge_predictions.R.
# ==============================================================================

RotationForest <- function(Df.X, Df.Y, kuse, juse, verbose = FALSE, ...) {
  RotForest <- list()
  class(RotForest) <- "RotationForest"
  fits <- list(); rots <- list()
  for (i in 1:juse) {
    model.current <- BuildOneModel(Df.X, Df.Y, kuse, ...)
    fits[[i]] <- model.current[[1]]
    rots[[i]] <- model.current[[2]]
    if (verbose) print(sprintf("Currently completed %i out of %i models", i, juse))
  }
  RotForest$models <- fits
  RotForest$rotations <- rots
  return(RotForest)
}

predict.RotationForest <- function(RotForest, Df.X, prob = FALSE) {
  prediction.probabilities <- list()
  for (i in 1:length(RotForest[[1]])) {
    model.current <- RotForest[[1]][[i]]
    data.current <- as.matrix(Df.X) %*% RotForest[[2]][[i]]
    data.current <- as.data.frame(data.current)
    colnames(data.current) <- paste0("X", 1:ncol(data.current))
    prediction.probabilities[[i]] <- as.matrix(predict(model.current, data.current))
  }
  results <- matrix(ncol = ncol(prediction.probabilities[[1]]), nrow = nrow(Df.X))
  colnames(results) <- colnames(prediction.probabilities[[1]])
  for (i in 1:nrow(Df.X)) {
    results[i, ] <- apply(do.call(rbind, lapply(prediction.probabilities, function(x) x[i, ])), 2, mean)
  }
  if (prob) return(results)
  return(apply(results, 1, function(x) names(which(x == max(x)))))
}

BuildOneModel <- function(Df.X, Df.Y, k, rows.use.frac = 0.75, ...) {
  M <- ceiling(ncol(Df.X) / k)
  R <- matrix(nrow = ncol(Df.X), ncol = ncol(Df.X), data = 0)
  R.order <- R
  Order <- data.frame(1:ncol(Df.X), sample(sort(rep(1:k, times = M))[1:ncol(Df.X)], size = ncol(Df.X), replace = FALSE))
  colnames(Order) <- c("V1", "V2")
  for (i in 1:k) {
    rows.use <- sample(1:nrow(Df.X), size = round(rows.use.frac * nrow(Df.X)), replace = FALSE)
    cols.use <- subset(Order, V2 == i)$V1
    start <- (i - 1) * M + 1
    end <- if (i != k) i * M else ncol(Df.X)
    Df.X.sub <- Df.X[rows.use, cols.use]
    Df.X.sub.rotation <- prcomp(Df.X.sub)$rotation
    R[start:end, start:end] <- Df.X.sub.rotation
    R.order[start:end, cols.use] <- R[start:end, start:end]
  }
  Df.X.rotate <- as.matrix(Df.X) %*% R.order
  Df.rotate.full <- data.frame(Df.Y, Df.X.rotate)
  colnames(Df.rotate.full)[1] <- "class"
  fit <- rpart::rpart(class ~ ., data = Df.rotate.full, ...)
  return(list(fit, R.order))
}

# Fits with k = ceiling(n_features/3), retrying with /5 if the first attempt
# errors (small groups can make prcomp() fail) -- same retry pattern used
# throughout the original benchmark scripts.
fit_rotation_forest <- function(x_train, y_train, ncol_ref) {
  den <- 3
  model <- tryCatch(
    RotationForest(x_train, y_train, kuse = max(1L, ceiling(ncol_ref / den)), juse = 100),
    error = function(e) { message("Failed with den=3, retrying with den=5..."); NULL }
  )
  if (is.null(model)) {
    den <- 5
    model <- RotationForest(x_train, y_train, kuse = max(1L, ceiling(ncol_ref / den)), juse = 100)
  }
  model
}
