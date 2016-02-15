\name{evaluateCorrection}
\alias{evaluateCorrection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quality assessment of a batch correction}
\description{Two criteria are employed to assess the quality of a batch
  correction. The first is based on a PCA decomposition; the other on
  the variation between (biological) replicates. This function is a
  wrapper for the two underlying functions - use these if you want to
  have access to more specific arguments.}
\usage{
evaluateCorrection(X, Y, what = c("duplo", "PCA"), ...)
}
\arguments{
  \item{X}{Data matrix, with samples in the rows and variables in the columns.}
  \item{Y}{A vector of batch labels.}
  \item{what}{Name of the quality criterion to be employed.}
  \item{\dots}{Further arguments to graphical functions.}
}
\value{A number indicating the relative success of the batch
  correction. For the PCA criterion, the number is the average
  Bhattacharya distance between the batches - smaller is better. For the
  duplo criterion, the number means the relative fraction of variation
  accounted for by the sample (in our example the genotype): bigger is
  better.
}
\author{Ron Wehrens}
\seealso{
  \code{\link{evaluateDuplos}}, \code{\link{evaluatePCA}}
}
\examples{
data(BC)
set.1.lod <- min(set.1[!is.na(set.1)])

## Figure 2 from the paper
par(mfrow = c(1,2))
## uncorrected data
huhnPCA <- evaluateCorrection(set.1, set.1.Y, what = "PCA",
                              plot = TRUE, legend.loc = "bottomright")
title(main = paste("Interbatch distance:", round(huhnPCA, 3)))

## do correction
set.1.corrected.Q <-
  apply(set.1, 2, doBC, ref.idx = which(set.1.Y$SCode == "ref"),
        batch.idx = set.1.Y$Batch, minBsamp = 4,
        seq.idx = set.1.Y$SeqNr, method = "lm",
        imputeVal = set.1.lod)
huhnPCA.A <- evaluateCorrection(set.1.corrected.Q, set.1.Y, what = "PCA",
                                plot = TRUE, legend.loc = "bottomright")
title(main = paste("Q: Interbatch distance:", round(huhnPCA.A, 3)))
}
\keyword{manip}
