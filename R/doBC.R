doBC <- function(Xvec, ref.idx, batch.idx, seq.idx,
                 result = c("correctedX", "corrections"),
                 method = c("lm", "rlm", "tobit"),
                 correctionFormula = formula("X ~ S * B"),
                 minBsamp = 4, imputeVal = NULL, ...) {
  result <- match.arg(result)
  method <- match.arg(method)

  if (is.null(imputeVal) & method == "tobit")
    stop("Tobit regression requires a value for 'imputeVal'")
  
  ## next line commented out to get rid of unused levels...
  ##   if (!is.factor(batch.idx))
  batch.idx <- factor(batch.idx)
  nbatches <- nlevels(batch.idx)
  if (is.factor(seq.idx)) seq.idx <- as.numeric(levels(seq.idx))[seq.idx]
  
  if (is.null(seq.idx) & method != "lm") {
    warning("Using method = 'lm' since seq.idx equals NULL")
  }

  ## convert TRUE/FALSE vector to vector of numerical indices
  if (is.logical(ref.idx)) ref.idx <- which(ref.idx)
  Xref <- Xvec[ref.idx]
  Bref <- batch.idx[ref.idx]
  Sref <- seq.idx[ref.idx]
  
  glMean <- mean(Xref, na.rm = TRUE)

  if (is.null(seq.idx)) {
    if (!is.null(imputeVal))
      Xref[is.na(Xref)] <- imputeVal
    
    Bmod <- lm(Xref ~ Bref - 1)
    Bcorrections <- (glMean - coef(Bmod))[ batch.idx ]

    switch(result,
           correctedX = Xvec + Bcorrections,
           Bcorrections)
  } else {
    ## check that at least minBsamp samples per batch are present
    ## if less than minBsamp samples present, remove batch for correction by
    ## setting values to NA. 
    nNonNA <- tapply(Xref, Bref, function(x) sum(!is.na(x)))
    tooFew <- names(nNonNA)[nNonNA < minBsamp]
    ## no batch found with enough samples...
    if (length(tooFew) > length(levels(Bref))-2)
      return(rep(NA, length(Xvec)))


    ## finally plug in imputeVal for non-detects
    if (!is.null(imputeVal))
      Xref[is.na(Xref)] <- imputeVal
    fitdf <- data.frame(S = Sref, B = Bref, X = Xref)
    ## subset argument for rlm does not work correctly: the number of
    ## levels is not diminished and as a result we get a singular
    ## matrix... the only solution is to take out the unused levels
    ## from the fitdf object.
    if (length(tooFew) > 0) {
      fitdf <- fitdf[!(fitdf$B %in% tooFew),]
      fitdf$B <- factor(fitdf$B)
    }
    
    Bmods2 <- switch(method,
                     lm = lm(correctionFormula, data = fitdf),
                     rlm = MASS:::rlm(correctionFormula, data = fitdf,
                                      maxit = 50),
                     tobit = AER:::tobit(correctionFormula, data = fitdf,
                             left = imputeVal))

    ## Now make predictions for each sample, not just the ref samples
    ## The actual correction is the following:
    ## ycorr = y - pred + gm
    predictdf <- data.frame(S = seq.idx, B = batch.idx)
    predictdf$B[predictdf$B %in% tooFew] <- NA
    predictions <- rep(NA, length(Xvec))
    predictions[!(predictdf$B %in% tooFew)] <-
      predict(Bmods2, newdata = predictdf)
 
    switch(result,
           correctedX = Xvec + glMean - predictions,
           glMean - predictions)
  }
}
