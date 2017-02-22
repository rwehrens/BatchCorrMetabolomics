evaluateCorrection <- function(X, Y, what = c("duplo", "PCA"), ...) {
  what <- match.arg(what)
  switch(what,
         duplo = evaluateDuplos(X, Y, ..., perMetabolite = FALSE),
         evaluatePCA(X, Y, ..., perBatch = FALSE))
}
