\name{BatchCorrMetabolomics-package}
\alias{BatchCorrMetabolomics-package}
\alias{BatchCorrMetabolomics}
\docType{package}
\title{
  Batch correlation for untargeted MS-based metabolomics data.
}
\description{Compendium to the paper
  "Improved batch correlation for untargeted MS-based metabolomics data"
  by Wehrens et al. (2016). All data are provided,
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
\references{
  "@Article{Wehrens2016,
  author =	 {Ron Wehrens and Jos.~A.~Hageman and Fred~van~Eeuwijk
                  and Rik~Kooke and P\'adraic~J.~Flood and Erik
                  Wijnker and Joost~J.B.~Keurentjes and Arjen~Lommen
                  and Henri\"ette~D.L.M.~van~Eekelen and
                  Robert~D.~Hall and Roland~Mumm and Ric~C.H.~de~Vos},
  title =	 {Improved batch correction in untargeted MS-based
                  metabolomics},
  journal =	 {Metabolomics},
  year =	 2016,
  volume =	 12,
  DOI =		 {10.1007/s11306-016-1015-8},
  pages =	 {1--12}

}"
}
\author{
  \packageAuthor{BatchCorrMetabolomics}
  
  Maintainer: \packageMaintainer{BatchCorrMetabolomics}
}
\keyword{ package }
\seealso{
  \code{\link{doBC}},  \code{\link{evaluateCorrection}}
}
