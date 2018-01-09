evaluatePCA <- function(X, Y, 
                        npc = 2, plot = FALSE, 
                        batch.colors, scaleX = TRUE,
                        legend.loc = "topright",
                        legend.col = 2, ..., perBatch = TRUE) {
  nbatches <- nlevels(Y$Batch)
  
  ## choose a default set of colors
  if (plot & missing(batch.colors)) {
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    batch.colors <- getPalette(nbatches)
  }

  noref.idx <- which(Y$SCode != "ref")
  Xsample <- X[noref.idx,]
  YSample <- Y[noref.idx,]
  
  Xsample <- Xsample[, apply(Xsample, 2, function(x) !all(is.na(x)))]
  
  ## replace NA values with column means
  for (i in 1:ncol(Xsample))
    Xsample[is.na(Xsample[,i]),i] <- mean(Xsample[,i], na.rm = TRUE)

  Xsample <- Xsample[, apply(Xsample, 2, sd, na.rm = TRUE) > 0]

  X.PCA <- PCA(scale(Xsample))
  
  if (plot) {
    scoreplot.PCA(X.PCA,
                  col = batch.colors[as.integer(Y$Batch)[noref.idx] ],
                  pch = as.integer(Y$Batch)[noref.idx], ...)
    legend(legend.loc, legend = levels(Y$Batch),
           col = batch.colors, pch = 1:nbatches, ncol = legend.col)
  }
  
  Xscores <- scores.PCA(X.PCA)[, 1:npc, drop = FALSE]
  ## a double loop is necessary to compare all batches... we
  ## only do the top triangle.
  batch.means <-
    lapply(levels(YSample$Batch),
           function(btch)
             colMeans(Xscores[which(YSample$Batch == btch),,drop=FALSE]))
  batch.covs <-
    lapply(levels(YSample$Batch),
           function(btch)
             cov(Xscores[which(YSample$Batch == btch),,drop=FALSE]))
  noCov.idx <- which(sapply(batch.covs, function(x) all(x < 1e-8)))
  if ((nnoCov <- length(noCov.idx)) > 0) {
    warning(paste("Too little information for batch correction in the following batches:\n",
                  levels(YSample$Batch)[noCov.idx],
                  "- ignoring these batches in the PCA criterion"))
    nbatches <- nbatches - nnoCov
    batch.covs <- batch.covs[-noCov.idx]
    batch.means <- batch.means[-noCov.idx]
  }
  
  batch.dist <- matrix(0, nbatches, nbatches)
  for (i in 2:nbatches)
    for (j in 1:(i-1))
      batch.dist[j, i] <- bhattacharyya.dist(batch.means[[j]],
                                             batch.means[[i]],
                                             batch.covs[[j]],
                                             batch.covs[[i]])

  if (perBatch) {
    batch.dist + t(batch.dist)
  } else {
    ## Here we take the mean and not the median since one deviating
    ## batch is already a problem.
    mean(batch.dist[col(batch.dist) > row(batch.dist)])
  }
}
