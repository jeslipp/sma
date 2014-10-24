### import libraries

library(limma)
library(impute)
library(GenomicFeatures)
library(dplyr)
library(broom)

### Main functions

#' Load raw data files using limma
#'
#' Loading, normalization, and generation of ExpressionSet using limma functions
#'
#' @param target_path character
#' @param norm_method character, default is "loess"
#' @param imputation character, default is "knn", other values will cause no imputation
#' @import limma
#' @import GenomicFeatures
#' @import impute
#' @return ExpressionSet, normalized Expression set
#' @export

sma.load <- function(target_path, norm_method = "loess", imputation = "knn") {
  # read target file
  targets <- readTargets(target_path)
  rownames(targets) <- targets$SampleName
  # read raw microarray data
  RG <- read.maimages(targets, 
                      source = "genepix", 
                      wt.fun = wtflags(weight=0, cutoff=-50), 
                      verbose = FALSE)
  # normalize raw microrray data
  MA <- normalizeWithinArrays(RG, method = norm_method)
  log2vals <- MA$M
  dimnames(log2vals) <- list(annotation$ArrayID, targets$SampleName)
  # impute missing values by k-nearest neighbour algorithm
  if (imputation == "knn") {
    log2vals <- impute.knn(log2vals)$data
  }
  # generate phenoData
  phenoData <- new("AnnotatedDataFrame", data = targets)
  # generate featureData
  rownames(annotation) <- annotation$ArrayID
  featureData <- new("AnnotatedDataFrame", data = annotation)
  # generate expression data
  experimentData <- new("MIAME",
                        name = "", 
                        lab = "", 
                        contact = "", 
                        title = "", 
                        abstract = "")
  
  # generate ExpressionSet
  es_raw <- new("ExpressionSet",
                exprs = log2vals,
                phenoData = phenoData, 
                experimentData = experimentData, 
                featureData = featureData, 
                annotation = "")
  return(es_raw)
}

#' Analyze splicing microarray data using E*I/J score
#'
#' Averaging of duplicate probes, calculation of score, and statistical analysis
#'
#' @param es ExpressionSet
#' @import limma
#' @import GenomicFeatures
#' @import dplyr
#' @return list, containing averaged ExpressionSet, score ExpressionSet, statistical analysis based in score ExpressionSet
#' @export

sma.analyze <- function(es) {
  ### build expression set with averaged probes
  log2_avg <- by(data = exprs(es), INDICES = fData(es)$ProbeName, FUN = colMeans)
  log2_avg <- do.call("rbind", log2_avg)
  idx <- match(rownames(log2_avg), fData(es)$ProbeName)
  featureData_avg <- fData(es)[idx, ]
  rownames(log2_avg) <- rownames(featureData_avg)
  featureData_avg <- new("AnnotatedDataFrame", data = featureData_avg)
  phenoData <- new("AnnotatedDataFrame", data = pData(es))
  es_avg <- new("ExpressionSet",
                exprs = log2_avg,
                phenoData = phenoData,
                experimentData = experimentData(es), 
                featureData = featureData_avg, 
                annotation = "")
  
  ### build expression set for score
  eij <- map_EIJ(fData(es_avg))
  eij <- eij[, 3:5]
  # calculate intron*exon/junction score
  log2_eij <- t(apply(eij, 1, calcScore, exprs(es_avg)))
  rownames(log2_eij) <- eij[, "Intron"]
  featureData_eij <- fData(es)[rownames(log2_eij), ]
  featureData_eij <- dplyr::select(featureData_eij, ArrayID, ProbeID, Gene, Systematic, 
                                   Position, ProbeName)
  featureData_eij <- new('AnnotatedDataFrame', data = featureData_eij)
  es_eij <- new("ExpressionSet",
                exprs = log2_eij,
                phenoData = phenoData,
                experimentData = experimentData(es),
                featureData = featureData_eij,
                annotation = "")
  
  ### calculate stats
  genotypes <- unique(pData(es_eij)$Genotype)
  stats <- lapply(genotypes, calc_stats, es_eij)
  names(stats) <- genotypes
  
  ### return results
  return(list(es_avg = es_avg, 
              es_eij = es_eij, 
              stats = stats))
}

#' Linear modeling of significantly retained introns
#'
#' Fits logistic regression model to classification into significantly retained and unaffected introns based on intron specific features
#'
#' @param stats list, containing statistical analysis obtained by sma.analyze
#' @param features data.frame, containing microarray specific intron features
#' @param p numeric, statistical cut-off to be used for analysis, default is FDR = 0.05
#' @import broom
#' @return list, containing logistic regression models
#' @export

sma.model <- function(stats, features, p = 0.05) {
  Y <- lapply(stats, classify_introns, p)
  Y <- lapply(Y, factor, levels = c("ns", "up", "down"), labels = c("ns", "up", "down"))
  names(Y) <- names(stats)
  Y <- Y[sapply(Y, function(y) sum(y == "up") > 500)]
  idx_match <- match(rownames(stats[[1]]), rownames(features))
  X <- features[idx_match, ]
  
  models <- list()
  for (genotype in names(Y)) {
    y <- Y[genotype]
    data <- cbind(y, X)
    colnames(data)[1] <- "class"
    data <- subset(data, class != "down" & complete.cases(data))
    mf <- glm(class ~ ., data = data, family = "binomial")
    cat(paste0("Simplifying model for ", genotype, "\n"))
    mr <- step(mf, trace = FALSE)
    models[[genotype]] <- tidy(mr)
  }
  return(models)
}

### Utility functions


# Map corresponding exons introns and junctions
map_EIJ <- function(featureData) {
  fd <- dplyr::select(featureData, ArrayID, Gene, Type, Position)
  ex <- fd %>% dplyr::select(-Position) %>% filter(Type == "Exon")
  it <- filter(fd, Type == "Intron")
  jc <- filter(fd, Type == "Junction")
  ij <- inner_join(it, jc, by = c("Gene", "Position"))
  eij <- inner_join(ij, ex, by = c("Gene"))
  eij <- dplyr::select(eij, Gene, Position, ArrayID, ArrayID.x, ArrayID.y)
  names(eij)[3:5] <- c("Exon", "Intron", "Junction")
  eij
}

# Calculate E*I/J score
calcScore <- function(row, log2vals, func) {
  apply(log2vals[row, ], 2, function(row) {row[1] + row[2] - row[3]})
}

# Statistical analysis using linear modeling followed by shrinkage of variance using empirical Bayes
calc_stats <- function(genotype, ExpressionSet, reference = "con") {
  genotypes <- pData(ExpressionSet)$Genotype
  idx <- which(genotypes == genotype)
  e <- ExpressionSet[, idx]
  dm <- modelMatrix(pData(e), ref = reference, verbose = FALSE)
  fit <- lmFit(exprs(e), dm)
  fit <- eBayes(fit)
  topTable(fit, adjust.method = "BH", number = Inf, sort.by = "none", genelist=fData(e)$PID)
}

classify_introns <- function(topTable, p = 0.05) {
  with(topTable, ifelse(adj.P.Val<=p & logFC>0, "up",
                        ifelse(adj.P.Val<=p & logFC<0, "down", "ns")))
}