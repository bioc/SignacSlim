## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning=FALSE, message=FALSE)

## -----------------------------------------------------------------------------
library(SignacSlim)

fpath <- system.file("extdata", "fragments.tsv.gz", package="SignacSlim")

ppath <- system.file("extdata", "peaks.rds", package="SignacSlim")

peaks <- readRDS(ppath)

fragments <- CreateFragmentObject(fpath)

Matrix <- FeatureMatrix(
  fragments = fragments,
  features = peaks
)

print(Matrix[1:10, 1:4])


## -----------------------------------------------------------------------------
sessionInfo()

