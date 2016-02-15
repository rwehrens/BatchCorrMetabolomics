evaluateDuplos <-
function(X, Y, 
                           plot = !perMetabolite, ## scaleX = TRUE,
                           perMetabolite = TRUE, ...) {
  nbatches <- nlevels(Y$Batch)
  
  Xsample <- X[Y$SCode != "ref",]
  Ysample <- Y[Y$SCode != "ref",]
  
  ## if (scaleX)
  ##   X <- scale(Xsample,
  ##              center = colMeans(Xsample, na.rm = TRUE),
  ##              scale = apply(Xsample, 2, sd, na.rm = TRUE))

  ## calcRep <- function(x, label) {
  ##   ivar <- aggregate(x, list(label), var, na.rm = TRUE)$x
  ##   count <- aggregate(x, list(label), length)$x - 1
  ##   if (all(is.na(ivar))) return(NA)
  ##   if (max(count) < 2) return(NA)
  ##   withinVar <- sum(ivar * count, na.rm = TRUE) / sum(count)
    
  ##   betweenVar <- var(aggregate(x, list(label), mean, na.rm = TRUE)$x)

  ##   betweenVar / (betweenVar + withinVar)
  ## }
  ## repeatability <- 
  ##   sapply(1:ncol(Xsample),
  ##          function(ii)
  ##          calcRep(Xsample[,ii], Ysample$SCode))


  repMats <- by(Xsample, list(factor(Ysample$SCode)), I)
  short.idx <- which(sapply(repMats, nrow) < 2)

  ## between-replicate variability
  cMeans <- sapply(repMats[-short.idx], colMeans, na.rm = TRUE)
  cMeans[!is.finite(cMeans)] <- NA
  bVars <- apply(cMeans, 1, var, na.rm = TRUE)

  ## within-replicates variability
  df <- sapply(repMats,
               function(x)
                 apply(x, 2, function(xx) sum(!is.na(xx)))) - 1
  metwVars <- sapply(repMats,
                     function(x)
                       apply(x, 2, var))
  wVars <- rowSums(metwVars * df, na.rm = TRUE) /
    apply(df, 1, function(x) sum(x[x > 0]))
  ## sometimes a metabolite is not present often enough in any group
  wVars[!is.finite(wVars)] <- NA

  repeatability <- bVars / (bVars + wVars)
  
  ## don't use the "main" argument in the function call, add
  ## the plot title later
  if (plot) {
    h2 <- hist(repeatability, main = "", xlab = "Repeatability", ...)
    abline(v = mean(repeatability, na.rm = TRUE), col = 2, lty = 2)
  }

  if (perMetabolite) {
    repeatability
  } else {
    mean(repeatability, na.rm = TRUE)
  }
}
