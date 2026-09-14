# ==============================================================================
# ODRF predict methods -- rewritten from scratch.
#
# The CRAN `ODRF` package's own `predict.ODRF`/predict.ODT did not work
# correctly for this benchmark's regression setting (split = "mse"), so these
# are hand-written replacements built directly on the package's internal tree
# structure (including its C entry point `_ODRF_predict_ODT`). Sourced by
# both run_benchmark_ODRF.R and run_benchmark_SPORF.R (SPORF = the same ODRF
# package with NodeRotateFun = "RotMatRand").
# ==============================================================================

predictTree <- function(structure, Xnew, split, Levels, ...) {
  nodeNumLabel <- structure$nodeNumLabel
  if (split %in% c("gini", "entropy")) {
    if (all(structure$nodeCutValue == 0))
      nodeNumLabel <- matrix(nodeNumLabel, nrow = 1, ncol = length(Levels))
    nodeLabel <- Levels[max.col(nodeNumLabel)]
    nodeLabel[which(rowSums(nodeNumLabel) == 0)] <- "0"
  } else {
    nodeLabel <- as.character(nodeNumLabel[, 1])
  }
  if (all(structure$nodeCutValue == 0)) {
    pred <- rep(nodeLabel, nrow(Xnew))
    node <- rep(0, nrow(Xnew))
  } else {
    predict_tree <- .Call("_ODRF_predict_ODT",
                          PACKAGE = "ODRF", Xnew, structure$nodeRotaMat,
                          structure$nodeCutValue, structure$childNode, nodeLabel)
    pred <- predict_tree$prediction
    node <- as.integer(predict_tree$node)
  }
  if (!split %in% c("gini", "entropy")) pred <- as.numeric(pred)
  return(list(prediction = pred, leafnode = node))
}

predict.ODT <- function(object, Xnew, leafnode = FALSE, ...) {
  pp <- object$data$p
  if (!is.null(object$data$catLabel) && (sum(object$data$Xcat) > 0))
    pp <- pp - length(unlist(object$data$catLabel)) + length(object$data$Xcat)
  if (ncol(Xnew) != pp) stop("The dimensions of 'Xnew' and training data do not match")
  Xna <- is.na(Xnew)
  if (any(Xna)) {
    xj <- which(colSums(Xna) > 0)
    for (j in xj) Xnew[Xna[, j], j] <- mean(Xnew[, j], na.rm = TRUE)
  }
  Xnew <- as.matrix(Xnew)
  p <- ncol(Xnew); n <- nrow(Xnew)
  Xcat <- object$data$Xcat; catLabel <- object$data$catLabel; numCat <- 0
  if (sum(Xcat) > 0) {
    xj <- 1; Xnew1 <- matrix(0, nrow = n, ncol = length(unlist(catLabel)))
    for (j in seq_along(Xcat)) {
      catMap <- which(catLabel[[j]] %in% unique(Xnew[, Xcat[j]]))
      indC <- catLabel[[j]][catMap]
      Xnewj <- (matrix(Xnew[, Xcat[j]], n, length(indC)) == matrix(indC, n, length(indC), byrow = TRUE)) + 0
      if (length(indC) > length(catLabel[[j]])) Xnewj <- Xnewj[, seq_along(catLabel[[j]])]
      xj1 <- xj + length(catLabel[[j]]); Xnew1[, (xj:(xj1 - 1))[catMap]] <- Xnewj; xj <- xj1
    }
    Xnew <- cbind(Xnew1, Xnew[, -Xcat]); p <- ncol(Xnew); numCat <- length(unlist(catLabel))
    rm(Xnew1, Xnewj)
  }
  if (!is.numeric(Xnew)) Xnew <- apply(Xnew, 2, as.numeric)
  if (object$data$Xscale != "No") {
    indp <- (numCat + 1):p
    Xnew[, indp] <- (Xnew[, indp] - matrix(object$data$minCol, n, length(indp), byrow = TRUE)) /
      matrix(object$data$maxminCol, n, length(indp), byrow = TRUE)
  }
  if (object$data$TreeRandRotate)
    Xnew[, object$data$rotdims] <- Xnew[, object$data$rotdims, drop = FALSE] %*% object$data$rotmat
  predict_tree <- predictTree(object$structure, Xnew, object$split, object$Levels)
  if (leafnode) predict_tree$leafnode else predict_tree$prediction
}

predict.ODRF <- function(object, Xnew, type = "response", weight.tree = FALSE, ...) {
  pp <- object$data$p
  if (!is.null(object$data$catLabel) && (sum(object$data$Xcat) > 0))
    pp <- pp - length(unlist(object$data$catLabel)) + length(object$data$Xcat)
  if (ncol(Xnew) != pp) stop("The dimensions of 'Xnew' and training data do not match")
  Xna <- is.na(Xnew)
  if (any(Xna)) {
    xj <- which(colSums(Xna) > 0)
    for (j in xj) Xnew[Xna[, j], j] <- mean(Xnew[, j], na.rm = TRUE)
  }
  Xnew <- as.matrix(Xnew); p <- ncol(Xnew); n <- nrow(Xnew)
  nC <- length(object$Levels); ntrees <- length(object$structure)
  Xcat <- object$data$Xcat; catLabel <- object$data$catLabel; numCat <- 0
  if (sum(Xcat) > 0) {
    xj <- 1; Xnew1 <- matrix(0, nrow = n, ncol = length(unlist(catLabel)))
    for (j in seq_along(Xcat)) {
      catMap <- which(catLabel[[j]] %in% unique(Xnew[, Xcat[j]]))
      indC <- catLabel[[j]][catMap]
      Xnewj <- (matrix(Xnew[, Xcat[j]], n, length(indC)) == matrix(indC, n, length(indC), byrow = TRUE)) + 0
      if (length(indC) > length(catLabel[[j]])) Xnewj <- Xnewj[, seq_along(catLabel[[j]])]
      xj1 <- xj + length(catLabel[[j]]); Xnew1[, (xj:(xj1 - 1))[catMap]] <- Xnewj; xj <- xj1
    }
    Xnew <- cbind(Xnew1, Xnew[, -Xcat]); p <- ncol(Xnew); numCat <- length(unlist(catLabel))
    rm(Xnew1, Xnewj)
  }
  if (!is.numeric(Xnew)) Xnew <- apply(Xnew, 2, as.numeric)
  if (object$data$Xscale != "No") {
    indp <- (sum(numCat) + 1):p
    Xnew[, indp] <- (Xnew[, indp] - matrix(object$data$minCol, n, length(indp), byrow = TRUE)) /
      matrix(object$data$maxminCol, n, length(indp), byrow = TRUE)
  }
  split <- object$split; Levels <- object$Levels; Rotate <- object$data$TreeRandRotate
  VALUE <- rep(ifelse(split == "mse", 0, "0"), n)
  TreePrediction <- vapply(object$structure, function(tree) {
    XXnew <- Xnew
    if (Rotate) XXnew[, tree$rotdims] <- XXnew[, tree$rotdims, drop = FALSE] %*% tree$rotmat
    predictTree(tree, XXnew, split, Levels)$prediction
  }, VALUE)
  Votes <- t(TreePrediction)
  weights <- rep(1, ntrees)
  if (weight.tree) {
    if (object$forest$ratOOB == 0) warning("ratOOB=0, weight.tree=TRUE invalid")
    else { oobErr <- sapply(object$structure, function(t) t$oobErr); weights <- 1 / (oobErr + 1e-5) }
  }
  weights <- weights / sum(weights)
  if (split != "mse") {
    weights <- rep(weights, n)
    Votes <- factor(c(Votes), levels = Levels)
    Votes <- as.integer(Votes) + nC * rep(0:(n - 1), rep(ntrees, n))
    Votes <- aggregate(c(rep(0, n * nC), weights), by = list(c(1:(n * nC), Votes)), sum)[, 2]
    prob <- matrix(Votes, n, nC, byrow = TRUE)
    prob <- prob / matrix(rowSums(prob), n, nC); colnames(prob) <- Levels
    pred <- Levels[max.col(prob)]
  } else {
    pred <- t(Votes) %*% weights
  }
  if (type == "response") return(pred)
  if (type == "prob")     return(prob)
  if (type == "tree")     return(TreePrediction)
}
