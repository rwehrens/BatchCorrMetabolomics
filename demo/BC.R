## Example R script reproducing all results in "Improved batch
## correction in untargeted MS-based metabolomics" by R. Wehrens et
## al., Metabolomics 12:1-12 (2016)
cat("Performing batch correction for three Arabidopsis metabolomics",
    "\ndata sets. In some cases, calculations may take a while - please",
    "\nbe patient...\n\n")

###################################################
### Load packages (quietly!) and example data
###################################################

suppressMessages(require(BatchCorrMetabolomics))
suppressMessages(require(RUVSeq))
data(BC)

###################################################
### data characteristics and definition of parameters
###################################################

minBatchOccurrence.Ave <- 2
minBatchOccurrence.Line <- 4

set.1.lod <- min(set.1[!is.na(set.1)])
set.2.lod <- min(set.2[!is.na(set.2)])
set.3.lod <- min(set.3[!is.na(set.3)])

LC.ngenotypes <- nlevels(set.1.Y$SCode)-1 ## take away the ref samples
LC.nref <- sum(set.1.Y$SCode == "ref")
LC.nNA <- apply(set.1, 2, function(x) sum(is.na(x)))

GC.ngenotypes <- nlevels(set.2.Y$SCode)-1 ## take away the ref samples
GC.nref <- sum(set.2.Y$SCode == "ref")
GC.nNA <- apply(set.2, 2, function(x) sum(is.na(x)))

TOF.ngenotypes <- nlevels(set.3.Y$SCode)-1 ## take away the ref samples
TOF.nNA <- apply(set.3, 2, function(x) sum(is.na(x)))
TOF.nref <- sum(set.3.Y$SCode == "ref")



###################################################
### Correction of set.1 - leads to object allResults
###################################################
cat("\nPerforming corrections for LC-MS data...")
conditions <- c("", "0", "1", "2", "c")
experiments <- c(t(outer(c("Q", "S"), conditions, paste, sep = "")))
methods <- rep("lm", length(experiments))
methods[grep("c", experiments)] <- "tobit"
imputeValues <- rep(NA, length(experiments))
imputeValues[grep("0", experiments)] <- 0
imputeValues[grep("1", experiments)] <- set.1.lod / 2
imputeValues[grep("2", experiments)] <- set.1.lod
imputeValues[grep("c", experiments)] <- set.1.lod - .01
refSamples <- list("Q" = which(set.1.Y$SCode == "ref"), 
                   "S" = which(set.1.Y$SCode != "ref"))
strategies <- rep(c("Q", "S"), each = length(conditions))

huhnDuplo <- evaluateCorrection(set.1, set.1.Y, what = "duplo", plot = FALSE)
huhnPCA <- evaluateCorrection(set.1, set.1.Y, what = "PCA", plot = FALSE)

## ANCOVA-type corrections
suppressMessages(allResults <-
                   lapply(seq(along = experiments),
                          function(ii)
                            apply(set.1, 2, doBC,
                                  ref.idx = refSamples[[ strategies[[ii]] ]],
                                  batch.idx = set.1.Y$Batch,
                                  minBsamp = minBatchOccurrence.Line,
                                  seq.idx = set.1.Y$SeqNr,
                                  method = methods[ii],
                                  imputeVal = imputeValues[ii])))
names(allResults) <- experiments

## add RUV results
idx <- which(set.1.Y$SCode == "ref")
replicates.ind <- matrix(-1, nrow(set.1) - length(idx) + 1, length(idx))
replicates.ind[1,] <- idx
replicates.ind[-1,1] <- (1:nrow(set.1))[-idx]

C0 <- C1 <- C2 <- set.1
C0[is.na(C0)] <- 0              #impute 0
C1[is.na(C1)] <- 0.5*set.1.lod  #impute 0.5*LOD
C2[is.na(C2)] <- set.1.lod      #impute LOD
nColumns <- ncol(C2)
w <- c('C0', 'C1', 'C2')

## RUVs gives a warning if a matrix contains negative numbers - it
## then assumes the matrix is already log-scaled. This is indeed the
## case here, so we suppress the warnings.
suppressWarnings(RUVresults <- 
                   lapply(w,
                          function(woppa) {
                            XimpRUV <- t(RUVs(t(get(woppa)),
                                              1:nColumns,
                                              k = 3,
                                              replicates.ind,
                                              round = FALSE)$normalizedCounts)
                            ## remove imputed values
                            XimpRUV[is.na(set.1)] <- NA
                            
                            log(XimpRUV)
                          }))
names(RUVresults) <- paste("R", 0:2, sep = "")
allResults <- c(allResults, RUVresults)


###################################################
### Correction of set.2 - objects GCMSResults.Ave, GCMSResults.Line,
### GCRUVResults
###################################################
cat("\nPerforming corrections for GC-QQQ-MS data...")
conditions <- c("", "0", "1", "2", "c")
experiments <- c(t(outer(c("Q", "S"), conditions, paste, sep = "")))
methods <- rep("lm", length(experiments))
methods[grep("c", experiments)] <- "tobit"
imputeValues <- rep(NA, length(experiments))
imputeValues[grep("0", experiments)] <- 0
imputeValues[grep("1", experiments)] <- set.2.lod / 2
imputeValues[grep("2", experiments)] <- set.2.lod
imputeValues[grep("c", experiments)] <- set.2.lod - .01
refSamples <- list("Q" = which(set.2.Y$SCode == "ref"),
                   "S" = which(set.2.Y$SCode != "ref"))
strategies <- rep(c("Q", "S"), each = length(conditions))

## leave out the censored regressions for A
exp.idx <- (1:length(experiments))[-grep("c", experiments)]
exp.idx.line <- which(strategies == "S")

suppressMessages(GCMSResults.Ave <-
                   lapply(exp.idx,
                          function(ii) 
                            apply(set.2, 2, doBC,
                                  ref.idx = refSamples[[ strategies[[ii]] ]],
                                  batch.idx = set.2.Y$Batch,
                                  minBsamp = minBatchOccurrence.Ave,
                                  seq.idx = set.2.Y$SeqNr,
                                  method = methods[ii],
                                  correctionFormula = formula("X ~ B"),
                                  imputeVal = imputeValues[ii])))
names(GCMSResults.Ave) <- experiments[exp.idx]

suppressMessages(GCMSResults.Line <-
                   lapply(exp.idx.line,
                          function(ii) 
                            apply(set.2, 2, doBC,
                                  ref.idx = refSamples[[ strategies[[ii]] ]],
                                  batch.idx = set.2.Y$Batch,
                                  minBsamp = minBatchOccurrence.Line,
                                  seq.idx = set.2.Y$SeqNr,
                                  method = methods[ii],
                                  imputeVal = imputeValues[ii])))
names(GCMSResults.Line) <- experiments[exp.idx.line]

## Application of RUV to these data:
idx <- which(set.2.Y$SCode == "ref")
replicates.ind <- matrix(-1, nrow(set.2) - length(idx) + 1, length(idx))
replicates.ind[1,] <- idx
replicates.ind[-1,1] <- (1:nrow(set.2))[-idx]

## prepare
C0 <- C1 <- C2 <- set.2
C0[is.na(C0)] <- 1e-5              # impute 0
C1[is.na(C1)] <- 0.5*set.2.lod     # impute 0.5*LOD
C2[is.na(C2)] <- set.2.lod         # impute LOD
nColumns <- ncol(C2)

w <- c('C0', 'C1', 'C2')
nw <- length(w)

suppressWarnings(GCRUVresults <-
                   lapply(w,
                          function(woppa) {
                            XimpRUV <- t(RUVs(t(get(woppa)),
                                              1:nColumns,
                                              k = 3,
                                              replicates.ind,
                                              round = FALSE)$normalizedCounts)
                            XimpRUV[is.na(set.2)] <- NA
                            log(XimpRUV)
                          }))

names(GCRUVresults) <- paste("R", 0:2, sep = "")

  
###################################################
### Correction of set.3 - objects allResultsAve, allResultsLine, allResultsRUV
###################################################
cat("\nPerforming corrections for GC-ToF-MS data...")
conditions <- c("", "0", "1", "2", "c")
experiments <- c(t(outer(c("Q", "S"), conditions, paste, sep = "")))
methods <- rep("lm", length(experiments))
methods[grep("c", experiments)] <- "tobit"
imputeValues <- rep(NA, length(experiments))
imputeValues[grep("0", experiments)] <- 0
imputeValues[grep("1", experiments)] <- set.3.lod / 2
imputeValues[grep("2", experiments)] <- set.3.lod
imputeValues[grep("c", experiments)] <- set.3.lod - .01
refSamples <- list("Q" = which(set.3.Y$SCode == "ref"),
                   "S" = which(set.3.Y$SCode != "ref"))
strategies <- rep(c("Q", "S"), each = length(conditions))

## leave out the censored regressions for A
exp.idx <- (1:length(experiments))[-grep("c", experiments)]
suppressMessages(allResultsAve <-
                   lapply(exp.idx,
                          function(ii)
                            apply(set.3, 2, doBC,
                                  ref.idx = refSamples[[ strategies[[ii]] ]],
                                  batch.idx = set.3.Y$Batch,
                                  minBsamp = minBatchOccurrence.Ave,
                                  correctionFormula = formula("X ~ B"),
                                  seq.idx = set.3.Y$SeqNr,
                                  method = methods[ii],
                                  imputeVal = imputeValues[ii])))
names(allResultsAve) <- experiments[exp.idx]

suppressMessages(allResultsLine <-
                   lapply(seq(along = experiments),
                          function(ii)
                            apply(set.3, 2, doBC,
                                  ref.idx = refSamples[[ strategies[[ii]] ]],
                                  batch.idx = set.3.Y$Batch,
                                  minBsamp = minBatchOccurrence.Line,
                                  seq.idx = set.3.Y$SeqNr,
                                  method = methods[ii],
                                  imputeVal = imputeValues[ii])))
names(allResultsLine) <- experiments

idx <- which(set.3.Y$SCode == "ref")
replicates.ind <- matrix(-1, nrow(set.3) - length(idx) + 1, length(idx))
replicates.ind[1,] <- idx
replicates.ind[-1,1] <- (1:nrow(set.3))[-idx]

suppressWarnings(allResultsRUV <-
                   lapply(0:2,
                          function(ImpVal) {
                            huhn <- set.3
                            huhn[is.na(huhn)] <- ImpVal * set.3.lod / 2
                            woppa <- t(RUVs(t(huhn),
                                            1:ncol(huhn),
                                            k = 3,
                                            replicates.ind,
                                            round = FALSE)$normalizedCounts)
                            woppa[is.na(set.3)] <- NA
                            log(woppa)
                          }))
  


###################################################
### Figure 1
###################################################
graphics.off()
dev.new(width = 14)
suppressMessages(example(doBC))
readline("Hit return to continue...")

###################################################
### Figure 2
###################################################
par(mfrow = c(1,2))
huhnPCA <- evaluateCorrection(set.1, set.1.Y, what = "PCA",
                              plot = TRUE, legend.loc = "bottomright")
title(main = paste("Interbatch distance:", round(huhnPCA, 3)))

huhnPCA.A <- evaluateCorrection(allResults[["Q"]], set.1.Y, what = "PCA",
                                plot = TRUE, legend.loc = "bottomright")
title(main = paste("Q: Interbatch distance:", round(huhnPCA.A, 3)))
readline("Hit return to continue...")

###################################################
### Do evaluation for all results
###################################################
## set.1 criteria - results
cat("\nPerforming evaluation of all corrections...")
results <- matrix(0, length(allResults) + 1, 2)
dimnames(results) <- list(c("No correction", names(allResults)),
                          c("PCA", "duplo"))
results["No correction",] <- c(huhnPCA, huhnDuplo)
for (exp in names(allResults)) {
  x <- set.1
  woppa <- allResults[[exp]]
  x[!is.na(woppa)] <- woppa[!is.na(woppa)]
  results[exp, 1] <- evaluateCorrection(x, set.1.Y, "PCA", plot = FALSE)
  results[exp, 2] <- evaluateCorrection(x, set.1.Y, "duplo", plot = FALSE)
}
huhnDuplo.A <- results["Q", 2]

## set.2 criteria - uncorrectedResults, gcmsresults.Ave,
## gcmsresults.Line, gcRUVresults
uncorrectedResults <-
  c(evaluateCorrection(set.2, set.2.Y, "PCA", plot = FALSE), 
    evaluateCorrection(set.2, set.2.Y, "duplo", plot = FALSE))

gcmsresults.Ave <-
  t(sapply(GCMSResults.Ave,
           function(x) {
             woppa <- set.2
             woppa[!is.na(x)] <- x[!is.na(x)]
             c(evaluateCorrection(woppa, set.2.Y, "PCA", plot = FALSE),
               evaluateCorrection(woppa, set.2.Y, "duplo", plot = FALSE))}))

gcmsresults.Line <-
  t(sapply(GCMSResults.Line,
           function(x) {
             woppa <- set.2
             woppa[!is.na(x)] <- x[!is.na(x)]
             c(evaluateCorrection(woppa, set.2.Y, "PCA", plot = FALSE),
               evaluateCorrection(woppa, set.2.Y, "duplo", plot = FALSE))}))

## woppa stuff not necessary for RUV
gcRUVresults <-
  t(sapply(GCRUVresults,
           function(x)
             c(evaluateCorrection(x, set.2.Y, what = "PCA",
                                  plot = FALSE),
               evaluateCorrection(x, set.2.Y, what = "duplo",
                                  plot = FALSE))))

## set.3 criteria - results.ave, results.line, results.ruv
results.ave <- results.line <- matrix(NA, length(experiments) + 1, 2)
dimnames(results.line) <- dimnames(results.ave) <-
  list(c("No correction", experiments), c("PCA", "duplo"))

results.line[1,1] <- results.ave[1,1] <-
  evaluateCorrection(set.3, set.3.Y, what = "PCA", plot = FALSE)
results.line[1,2] <- results.ave[1,2] <-
  evaluateCorrection(set.3, set.3.Y, what = "duplo", plot = FALSE)

for (exp in experiments[exp.idx]) {
  x <- allResultsAve[[exp]]
  woppa <- set.3
  woppa[!is.na(x)] <- x[!is.na(x)]
  results.ave[exp, 1] <-
    evaluateCorrection(woppa, set.3.Y, "PCA", plot = FALSE)
  results.ave[exp, 2] <-
    evaluateCorrection(woppa, set.3.Y, "duplo", plot = FALSE)
}

for (exp in experiments) {
  x <- allResultsLine[[exp]]
  woppa <- set.3
  woppa[!is.na(x)] <- x[!is.na(x)]
  results.line[exp, 1] <-
    evaluateCorrection(woppa, set.3.Y, "PCA", plot = FALSE)
  results.line[exp, 2] <-
    evaluateCorrection(woppa, set.3.Y, "duplo", plot = FALSE)
}

results.ruv <-
  cbind(sapply(allResultsRUV,
               evaluateCorrection, set.3.Y,
               what = "PCA", plot = FALSE),
        sapply(allResultsRUV,
               evaluateCorrection, set.3.Y,
               what = "duplo", plot = FALSE))
dimnames(results.ruv) <- list(paste("R", 0:2, sep = ""),
                              c("PCA", "duplo"))


###################################################
### Figure 3
###################################################
huhnDuplo.metab <- evaluateDuplos(set.1, set.1.Y, plot = FALSE, 
                                  perMetabolite = TRUE)
woppa <- set.1
woppa[!is.na(allResults[["Q"]])] <- allResults[["Q"]][!is.na(allResults[["Q"]])]
      
huhnDuplo.metabA <- evaluateDuplos(woppa, set.1.Y, plot = FALSE, 
                                   perMetabolite = TRUE)
xyranges <- range(c(huhnDuplo.metab, huhnDuplo.metabA), na.rm = TRUE)
plot(huhnDuplo.metab, huhnDuplo.metabA, 
     main = "Metabolite repeatabilities",
     xlim = xyranges, ylim=xyranges,
     xlab = "Uncorrected data", ylab = "Corrected data (Q)")
abline(0, 1, col = "gray")
readline("Hit return to continue...")


###################################################
### Figure 4
###################################################
results.label <- factor(substr(rownames(results), 1, 1))
xlim.small <- c(.06, .135)
ylim.small <- c(.61, .65)
par(mfrow = c(1,2))
plot(results, main = "Data set I", 
     xlab = "Interbatch distance", ylab = "Repeatability",
     col = as.integer(results.label))
text(results, labels = rownames(results), 
     pos = ifelse(results.label == "N", 2, 4),
     col = as.integer(results.label))
plot(results, main = "Data set I (zoom)", xlim = xlim.small,
     xlab = "Interbatch distance", ylab = "Repeatability",
     ylim = ylim.small, col = as.integer(results.label))
text(results, labels = rownames(results), 
     pos = ifelse(results.label == "R", 2, 4),
     col = as.integer(results.label))
readline("Hit return to continue...")



###################################################
### Figure 5
###################################################
gcresults.ave <- rbind(uncorrectedResults,
                       gcmsresults.Ave,
                       gcRUVresults)
rownames(gcresults.ave)[1] <- "No correction"
colnames(gcresults.ave) <- c("PCA", "duplo")
gcresults.ave.label <- factor(substr(rownames(gcresults.ave), 1, 1))

gcresults.line <- rbind(uncorrectedResults,
                        gcmsresults.Line,
                        gcRUVresults)
rownames(gcresults.line)[1] <- "No correction"
colnames(gcresults.line) <- c("PCA", "duplo")

xl <- range(c(gcresults.line[,1], gcresults.ave[,1]))
yl <- range(c(gcresults.line[,2], gcresults.ave[,2]))

par(mfrow = c(1,2))
plot(gcresults.ave[,1], gcresults.ave[,2], 
     main = "Data set II - only batch correction",
     ylab= "Repeatability", xlab = "Interbatch distance",
     col = as.integer(gcresults.ave.label),
     ylim = yl, xlim = xl)
text(gcresults.ave[,1], gcresults.ave[,2], 
     labels = rownames(gcresults.ave), 
     pos = ifelse(gcresults.ave.label %in% c("N"), 2, 4),
     col = as.integer(gcresults.ave.label))

gcresults.line.label <- factor(substr(rownames(gcresults.line), 1, 1),
                               levels = levels(gcresults.ave.label))
plot(gcresults.line[,1], gcresults.line[,2], 
     main = "Data set II - batch and order correction", 
     ylab= "Repeatability", xlab = "Interbatch distance",
     col = as.integer(gcresults.line.label),
     xlim = xl, ylim = yl)
text(gcresults.line[,1], gcresults.line[,2], 
     labels = rownames(gcresults.line), 
     pos = ifelse(gcresults.line.label %in% c("N"), 2, 4),
     col = as.integer(gcresults.line.label))
readline("Hit return to continue...")


###################################################
### Figure 6
###################################################
floodresults.ave <- rbind(results.ave, results.ruv)
floodresults.ave.label <- factor(substr(rownames(floodresults.ave), 1, 1))
floodresults.line <- rbind(results.line, results.ruv)
floodresults.line.label <- factor(substr(rownames(floodresults.line), 1, 1))
PCA.range <- range(c(floodresults.ave[,1], floodresults.line[,1]), 
                   na.rm = TRUE)
duplo.range <- range(c(floodresults.ave[,2], floodresults.line[,2]),
                     na.rm = TRUE)
par(mfrow = c(1,2))
plot(floodresults.ave[,1], floodresults.ave[,2],
     main = "Data set III - only batch correction",
     ylab= "Repeatability", xlab = "Interbatch distance",
     xlim = PCA.range, ylim = duplo.range, 
     col = as.integer(floodresults.ave.label))
text(floodresults.ave[,1], floodresults.ave[,2],
     pos = ifelse(floodresults.ave.label == "N", 2, 4),
     col = as.integer(floodresults.ave.label),
     labels = rownames(floodresults.ave))

plot(floodresults.line[,1], floodresults.line[,2], 
     main = "Data set III - batch and order correction", 
     xlim = PCA.range, ylim = duplo.range, 
     ylab= "Repeatability", xlab = "Interbatch distance",
     col = as.integer(floodresults.line.label))
text(floodresults.line[,1], floodresults.line[,2],
     pos = ifelse(floodresults.line.label %in% c("N", "Q"), 2, 4),
     col = as.integer(floodresults.line.label),
     labels = rownames(floodresults.line))
readline("Hit return to continue...")

cat("\nCreating Table 2 from the paper...")

###################################################
### Table 2
###################################################
batch.presence.raw <- 
  apply(set.1, 2, function(y) 
    aggregate(y, list(set.1.Y$Batch), function(x) sum(!is.na(x)))$x)
batch.presence.A <- 
  apply(allResults[["Q"]], 2, function(y) 
    aggregate(y, list(set.1.Y$Batch), function(x) sum(!is.na(x)))$x)
batch.presence.B <- 
  apply(allResults[["S"]], 2, function(y) 
      aggregate(y, list(set.1.Y$Batch), function(x) sum(!is.na(x)))$x)

diffA <- sum(batch.presence.raw - batch.presence.A > 0)
diffB <- sum(batch.presence.raw - batch.presence.B > 0)
totalCombI <- nlevels(set.1.Y$Batch)*ncol(set.1)


batch.presence.rawII <- 
  apply(set.2, 2, function(y) 
    aggregate(y, list(set.2.Y$Batch), function(x) sum(!is.na(x)))$x)
batch.presence.AII.ave <- 
  apply(GCMSResults.Ave[["Q"]], 2, function(y) 
    aggregate(y, list(set.2.Y$Batch), function(x) sum(!is.na(x)))$x)
batch.presence.BII.ave <- 
  apply(GCMSResults.Ave[["S"]], 2, function(y) 
    aggregate(y, list(set.2.Y$Batch), function(x) sum(!is.na(x)))$x)
batch.presence.BII.line <- 
  apply(GCMSResults.Line[["S"]], 2, function(y) 
    aggregate(y, list(set.2.Y$Batch), function(x) sum(!is.na(x)))$x)

diffAII.ave <- sum(batch.presence.rawII - batch.presence.AII.ave > 0)
diffBII.ave <- sum(batch.presence.rawII - batch.presence.BII.ave > 0)
diffBII.line <- sum(batch.presence.rawII - batch.presence.BII.line > 0)
totalCombII <- ncol(set.2) * nlevels(set.2.Y$Batch)


batch.presence.rawIII <- 
  apply(set.3, 2, function(y) 
    aggregate(y, list(set.3.Y$Batch), function(x) sum(!is.na(x)))$x)
batch.presence.AIII.ave <- 
  apply(allResultsAve[["Q"]], 2, function(y) 
    aggregate(y, list(set.3.Y$Batch), function(x) sum(!is.na(x)))$x)
batch.presence.AIII.line <- 
  apply(allResultsLine[["Q"]], 2, function(y) 
    aggregate(y, list(set.3.Y$Batch), function(x) sum(!is.na(x)))$x)
batch.presence.BIII.ave <- 
  apply(allResultsAve[["S"]], 2, function(y) 
    aggregate(y, list(set.3.Y$Batch), function(x) sum(!is.na(x)))$x)
batch.presence.BIII.line <- 
  apply(allResultsLine[["S"]], 2, function(y) 
    aggregate(y, list(set.3.Y$Batch), function(x) sum(!is.na(x)))$x)

diffAIII.ave <- sum(batch.presence.rawIII - batch.presence.AIII.ave > 0)
diffAIII.line <- sum(batch.presence.rawIII - batch.presence.AIII.line > 0)
diffBIII.ave <- sum(batch.presence.rawIII - batch.presence.BIII.ave > 0)
diffBIII.line <- sum(batch.presence.rawIII - batch.presence.BIII.line > 0)
totalCombIII <- ncol(set.3) * nlevels(set.3.Y$Batch)


countTab <- matrix(NA, 5, 3,
                   dimnames = list(c("Q (ave)", "Q (lin)", "S (ave)",
                                     "S (lin)", "R"),
                                   paste("Data set", c("I", "II", "III"))))
countTab[1, 2:3] <- c(100*diffAII.ave / totalCombII,
                      100*diffAIII.ave / totalCombIII)
countTab[2, c(1,3)] <- c(100*diffA / totalCombI,
                         100*diffAIII.line / totalCombIII)
countTab[3, 2:3] <- c(100*diffBII.ave / totalCombII,
                      100*diffBIII.ave / totalCombIII)
countTab[4,] <- c(100*diffB / totalCombI,
                  100*diffBII.line / totalCombII,
                  100*diffBIII.line / totalCombIII)
countTab[5,] <- 0

print(countTab, digits = 3)
