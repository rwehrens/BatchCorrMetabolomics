\name{evaluateDuplos}
\alias{evaluateDuplos}
\title{Evaluation of differences between duplicate samples}
\description{After batch correction, duplicate samples (either
  biological replicates or technical replicates) should be closer than
  in the uncorrected data. This function calculates the fraction of
  variance explained by the sample identity and the total variance, a
  measure also known as repeatability (or in genetic studies: heritability).
}
\usage{
evaluateDuplos(X, Y, plot = !perMetabolite, perMetabolite = TRUE, ...)
}
\arguments{
  \item{X}{Data matrix.}
  \item{Y}{Data.frame describing information of the samples: columns
    SCode and Batch are mandatory.}
  \item{plot}{Logical: should a histogram of repeatabilities be shown?}
  \item{perMetabolite}{Separate repeatabilities for each metabolite
    (variable), or one average?}
  \item{\dots}{Further graphical arguments.}
}
\value{The ratio of the amount of variation explained by the sample type
(here: genotype) and the total variation. This is also called
repeatability, or in a genetic context, heritability. Either given as
one overall average, or per sample type.}
\author{Ron Wehrens}

\seealso{
\code{\link{evaluatePCA}}
}
\examples{
data(BC)
set.1.lod <- min(set.1[!is.na(set.1)])

## do correction for the first 10 metabolites
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

huhnDuplos.A0 <- evaluateDuplos(set.1.corrected.Q0, set.1.Y,
                             plot = TRUE, perMetabolite = FALSE)
title(main = paste("Q: repeatability:", round(huhnDuplos.A0, 3)),
      sub = "NA imputation: 0")
huhnDuplos.A2 <- evaluateDuplos(set.1.corrected.Q2, set.1.Y,
                             plot = TRUE, perMetabolite = FALSE)
title(main = paste("Q: repeatability:", round(huhnDuplos.A2, 3)),
      sub = "NA imputation: LOD")

## Per-metabolite information:
evaluateDuplos(set.1.corrected.Q2, set.1.Y, what = "Duplos",
               perMetabolite = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip }
