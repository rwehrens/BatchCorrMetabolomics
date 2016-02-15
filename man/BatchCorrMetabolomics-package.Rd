\name{BatchCorrMetabolomics-package}
\alias{BatchCorrMetabolomics-package}
\alias{BatchCorrMetabolomics}
\docType{package}
\title{
  Batch correlation for untargeted MS-based metabolomics data.
}
\description{Compendium to the paper
  "Improved batch correlation for untargeted MS-based metabolomics data"
  by Wehrens et al., submitted for publication. All data are provided,
  as well as functions to reproduce the results in that paper. A demo R
  scripts is included that does this: simply type
  \code{demo("BC", echo = FALSE)}
  at the R prompt after loading the package.}
\details{
  The basic function in the package is \code{doBC} which does batch
  correction for one variable, given batch and measurement order
  information. Evaluation of batch correction for a number of variables
  simultaneously is available through \code{evaluateCorrection}.
}
\author{
  \packageAuthor{BatchCorrMetabolomics}
  
  Maintainer: \packageMaintainer{BatchCorrMetabolomics}
}
\keyword{ package }
\seealso{
  \code{\link{doBC}},  \code{\link{evaluateCorrection}}
}
