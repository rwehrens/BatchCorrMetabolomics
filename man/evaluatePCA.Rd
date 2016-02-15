\name{evaluatePCA}
\alias{evaluatePCA}
\title{Evaluation of differences between batches, using PCA}
\description{This function calculates the (average) Bhattacharyya
  distance between the batches of a data set. The lower this distance
  is, the more alike the batches. Alternatively, a distance matrix can
  be returned indicating for any pair of batches their B. distance.
}
\usage{
evaluatePCA(X, Y, npc = 2, plot = FALSE, batch.colors, scaleX = TRUE,
            legend.loc = "topright", legend.col = 2, ..., perBatch = TRUE)
}
\arguments{
  \item{X}{
    Data matrix: rows are samples, columns are features (metabolites in
    this case).}
  \item{Y}{Batch information: a data.frame with columns SCode, Batch and
  }
  \item{npc}{Number of PCs to include in the low-dimensional representation.}
  \item{plot}{Logical: should a score plot be shown?}
  \item{batch.colors}{Colors to be used for individual batches.}
  \item{scaleX}{Logical: should standardization (zero mean, unit
    variance) be applied for all columns? Default: yes.}
  \item{legend.loc}{Location of the legend.}
  \item{legend.col}{Number of columns in the legend.}
  \item{\dots}{Further graphical arguments.}
  \item{perBatch}{Logical: should the result be given as a distance
    matrix between batches (the default), or as one average distance?}
}
\value{Returns Bhattacharyya distances between batches. If
  \code{perBatch == TRUE}, a distance matrix is returned, otherwise the
  average value of all distances is returned.}
\author{Ron Wehrens}

\seealso{
\code{\link{evaluateDuplos}}
}

\examples{
data(BC)
set.1.lod <- min(set.1[!is.na(set.1)])

## do correction, only first ten metabolites of set.1
set.1.corrected.Q0 <-
  apply(set.1[,1:10], 2, doBC, ref.idx = which(set.1.Y$SCode == "ref"),
        batch.idx = set.1.Y$Batch, minBsamp = 4,
        seq.idx = set.1.Y$SeqNr, method = "lm",
        imputeVal = 0)
set.1.corrected.Q2 <-
  apply(set.1[,1:10], 2, doBC, ref.idx = which(set.1.Y$SCode == "ref"),
        batch.idx = set.1.Y$Batch, minBsamp = 4,
        seq.idx = set.1.Y$SeqNr, method = "lm",
        imputeVal = set.1.lod)

huhnPCA.A0 <- evaluatePCA(set.1.corrected.Q0, set.1.Y, perBatch = FALSE,
                          plot = TRUE, legend.loc = "bottomright")
title(main = paste("Q: Interbatch distance:", round(huhnPCA.A0, 3)),
      sub = "NA imputation: 0")
huhnPCA.A2 <- evaluatePCA(set.1.corrected.Q2, set.1.Y, perBatch = FALSE,
                          plot = TRUE, legend.loc = "bottomright")
title(main = paste("Q: Interbatch distance:", round(huhnPCA.A2, 3)),
      sub = "NA imputation: LOD")

## which batches are more similar?
B2B <- evaluatePCA(set.1.corrected.Q2, set.1.Y, what = "PCA", plot = FALSE,
                   perBatch = TRUE)
dimnames(B2B) <- list(levels(set.1.Y$Batch), levels(set.1.Y$Batch))
plot(hclust(as.dist(B2B)))
}
\keyword{manip}
