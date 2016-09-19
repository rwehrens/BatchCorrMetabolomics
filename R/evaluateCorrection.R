evaluateCorrection <- function(X, Y, what = c("duplo", "PCA"), ...) {
  what <- match.arg(what)
  switch(what,
         duplo = evaluateDuplos(X, Y, ..., plot = FALSE, perMetabolite = FALSE),
         evaluatePCA(X, Y, ..., plot = FALSE, perBatch = FALSE))
}
