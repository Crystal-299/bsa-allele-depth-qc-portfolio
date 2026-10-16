#!/usr/bin/env Rscript
# Optional dense-SNP analysis; requires QTLseqr. Not part of the tiny fixture test.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) stop("usage: Rscript scripts/qtlseqr_optional.R sites.tsv F2 high_size low_size window_bp outdir")
if (!requireNamespace("QTLseqr", quietly = TRUE)) stop("install QTLseqr from its official repository first")
input <- args[1]; population <- args[2]
high_size <- as.integer(args[3]); low_size <- as.integer(args[4]); window <- as.integer(args[5]); outdir <- args[6]
if (population != "F2") stop("This optional simulation is configured only for a verified F2 design")
if (any(is.na(c(high_size, low_size, window))) || min(high_size, low_size, window) < 1) stop("sizes/window must be positive integers")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
tab <- read.delim(input, check.names = FALSE)
if (!all(c("CHROM", "POS", "AD_REF.High", "AD_ALT.High", "AD_REF.Low", "AD_ALT.Low") %in% names(tab))) stop("missing QTLseqr columns")
chroms <- unique(tab$CHROM)
df <- QTLseqr::importFromTable(file = input, highBulk = "High", lowBulk = "Low", chromList = chroms)
df <- QTLseqr::filterSNPs(SNPset = df, refAlleleFreq = 0.15, minTotalDepth = 10,
                           maxTotalDepth = 500, minSampleDepth = 10)
if (nrow(df) < 20) stop("too few filtered SNPs for a meaningful windowed QTLseqr analysis")
df <- QTLseqr::runGprimeAnalysis(SNPset = df, windowSize = window, outlierFilter = "deltaSNP")
df <- QTLseqr::runQTLseqAnalysis(SNPset = df, windowSize = window, popStruc = "F2",
                                  bulkSize = c(high_size, low_size), replications = 1000,
                                  intervals = c(95, 99))
write.csv(df, file.path(outdir, "qtlseqr_sites.csv"), row.names = FALSE)
pdf(file.path(outdir, "qtlseqr_plots.pdf"), width = 11, height = 7)
print(QTLseqr::plotQTLStats(SNPset = df, var = "Gprime", plotThreshold = TRUE, q = 0.01))
print(QTLseqr::plotQTLStats(SNPset = df, var = "deltaSNP", plotIntervals = TRUE))
dev.off()
