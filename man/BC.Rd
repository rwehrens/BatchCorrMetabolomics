\name{BC}
\alias{BC}
\alias{set.1}
\alias{set.1.Y}
\alias{set.2}
\alias{set.2.Y}
\alias{set.3}
\alias{set.3.Y}
\docType{data}
\title{Three metabolomics data sets for Arabidopsis}
\description{Three real-life data sets are provided for testing batch
  correction procedures, all three from Arabidopsis thaliana samples. In
  each case, the data object (set.1, set.2 or set.3) contains relative
  intensities of reconstructed metabolites (columns) for a number of
  samples (rows). The corresponding meta-information object (set.1.Y,
  set.2.Y, set.3.Y, respectively) contains for every sample batch and
  sequence information, as well as (coded) information on the genotype
  and the sample type (study sample or reference sample). Processing of
  the raw data has been done with Metalign and MSClust programs.
}
\usage{data(BC)}
\details{Set 1: data from a set of 761 injections containing 51
  QC samples and biological replicates of 357 natural Arabidopsis
  thaliana accessions. Number of metabolites: 567 reconstructed
  metabolites.

  Set 2: data from a set of 753 injections containing 50
  QC samples and biological replicates of 357 natural Arabidopsis
  thaliana accessions (same as set.1). Number of metabolites: 603.

  Set 3: Data from a set of 240 injections containing 14
  QC samples and biological replicates of 58 Arabidopsis
  thaliana accessions. Number of metabolites: 75.
}
\references{
  @Article{Wehrens2016,
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

  }
\examples{
data(BC)
dim(set.1)
sum(is.na(set.1))
dim(set.1.Y)

dim(set.2)
sum(is.na(set.2))
dim(set.2.Y)

dim(set.3)
sum(is.na(set.3))
dim(set.3.Y)
}
\keyword{datasets}
