\name{doBC}
\alias{doBC}
\title{Batch correction for untargeted metabolomics data}
\description{
Variable-wise batch correction using explicit information on batch and
injection order. This is done through fitting ancova models with batch
as factorial information, and injection information as
numerical. Several types of regression are supported.
}
\usage{
doBC(Xvec, ref.idx, batch.idx, seq.idx,
     result = c("correctedX", "corrections"),
     method = c("lm", "rlm", "tobit"),
     correctionFormula = formula("X ~ S * B"),
     minBsamp = ifelse(is.null(seq.idx), 2, 4),
     imputeVal = NULL, ...)
}

\arguments{
  \item{Xvec}{A vector of intensity values (typically from one metabolite, but
    possibly also from one particular mass peak or (for NMR) chemical shift).}
  \item{ref.idx}{Indices for reference samples, i.e., samples that are
    expected to be stable during the measurement sequence.}
  \item{batch.idx}{Batch indices for all samples.}
  \item{seq.idx}{Information on within-batch measurement order: simply a
    vector of integers.}
  \item{result}{Determines whether the corrected data are returned
    (default), or the corrections themselves.}
  \item{method}{Defines the type of regression: default is to use simple
    least-squares regression ("lm"), but also robust regression ("rlm")
    and censored regression ("tobit") are supported.}
  \item{correctionFormula}{Model to be fitted in the regression step.}
  \item{minBsamp}{Minimum number of reference samples in each batch; if
    a batch has fewer reference samples, no correction is performed and a
    vector of NAs is returned.}
  \item{imputeVal}{What to impute for an NA value. If ImputeVal equals
    NULL (the default), NAs are retained.}
  \item{\dots}{additional plotting arguments.}
}
\value{The function returns a vector of the same size as the input
  data. The vector corresponds to either the corrected data, or to the
  corrections.
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
\author{Ron Wehrens}

\examples{
data(BC) ## provides set.1, set.2 and set.3

## Correct the second metabolite from set.1
metabo.idx <- 2
ref.idx <- which(set.1.Y$SCode == "ref")
M2c <- doBC(set.1[,metabo.idx], ref.idx = ref.idx,
            batch.idx = set.1.Y$Batch, method = "rlm",
            seq.idx = set.1.Y$SeqNr)
M2c <- doBC(set.1[,metabo.idx], ref.idx = ref.idx,
            batch.idx = set.1.Y$Batch, method = "tobit",
            seq.idx = set.1.Y$SeqNr, imputeVal = 0)
M2c <- doBC(set.1[,metabo.idx], ref.idx = ref.idx,
            batch.idx = set.1.Y$Batch,
            seq.idx = set.1.Y$SeqNr)

require(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
batch.colors <- getPalette(nlevels(set.1.Y$Batch))
plot(set.1[, metabo.idx], col = batch.colors[set.1.Y$Batch],
     main = "Raw values for metabolite 2 (set.1)")
plot(M2c, pch = 2, col = batch.colors[set.1.Y$Batch],
     main = "Corrected values for metabolite 2 (set.1)")

## Figure 1 from the paper: correction of metabolite 2 from set.1
par(mfrow = c(1,2))
batches <- c("B1", "B2")
sample.idx <- which(set.1.Y$Batch \%in\% batches)

plot(set.1.Y$SeqNr[sample.idx], set.1[sample.idx, metabo.idx], 
     col = as.numeric(set.1.Y[sample.idx, "SCode"] == "ref") + 1,
     pch = c(1, 19)[as.numeric(set.1.Y[sample.idx, "SCode"] == "ref") + 1],
     xlab = "Injection number", ylab = "Intensity (log-scaled)",
     main = paste("Metabolite", metabo.idx, "- before correction"))

## we assume consecutive batches
batch.lims <- aggregate(set.1.Y$SeqNr[sample.idx], 
                        by = list(set.1.Y$Batch[sample.idx]),
                        FUN = range)$x
abline(v = mean(batch.lims[1,2], batch.lims[2,1]), lty = 2)
text(aggregate(set.1.Y$SeqNr[sample.idx], 
               by = list(set.1.Y$Batch[sample.idx]), 
               FUN = mean)$x,
     par("usr")[4],
     batches, pos = 1, col = 4)

pred.df <- data.frame(y = set.1[sample.idx[ref.idx], metabo.idx],
                      x = set.1.Y$SeqNr[sample.idx[ref.idx]],
                      batch = set.1.Y$Batch[sample.idx[ref.idx]])
cormods <- lapply(batches,
                  function(b.idx) lm(y ~ x, subset = batch == b.idx,
                                     data = pred.df))
for (i in seq(along = cormods))
  lines(batch.lims[i,],
        predict(cormods[[i]], new = data.frame(x = batch.lims[i,])),
        col = 2)

## get corrected values
corvals <- doBC(set.1[sample.idx, metabo.idx],
                ref.idx = ref.idx,
                batch.idx = set.1.Y$Batch[sample.idx],
                seq.idx = set.1.Y$SeqNr[sample.idx], plot = FALSE)
plot(set.1.Y$SeqNr[sample.idx], corvals,
     col = as.numeric(set.1.Y[sample.idx, "SCode"] == "ref") + 1,
     pch = c(1, 19)[as.numeric(set.1.Y[sample.idx, "SCode"] == "ref") + 1],
     xlab = "Injection number", ylab = "Intensity (log-scaled)",
     main = paste("Metabolite", metabo.idx, "- after correction"))
abline(v = mean(batch.lims[1,2], batch.lims[2,1]), lty = 2)
text(aggregate(set.1.Y$SeqNr[sample.idx], by = list(set.1.Y$Batch[sample.idx]), 
               FUN = mean)$x,
     par("usr")[4],
     batches, pos = 1, col = 4)

abline(h = mean(corvals[ref.idx], na.rm = TRUE), col = 2)
}
\keyword{manip}
