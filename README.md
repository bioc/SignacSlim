# SignacSlim
A slim version of Signac for developing scATAC tools using signac functions.

## Why SignacSlim is developed?

Single cell ATAC-seq is a promising tool for analyzing transcriptional regulation. In addition, Seurat and Signac provides popular functions for single cell data pre-processing. However, Seurat and Signac depends on too many packages which makes them difficult for other package re-use. Signacslim is a slim version of Signac and it contains many useful methods from Signac and Seurat. SignacSlim depends less packages and can be easily installed copared with Signac and Seurat.


## Demo for using SignacSlim

Gnerating feature matrix for >10K single cells is time-consuming in R. Using "FeatureMatrix" function provided by SignacSlim will be much faster.

```
fpath <- system.file("extdata", "fragments.tsv.gz", package="SignacSlim")

ppath <- system.file("extdata", "peaks.rds", package="SignacSlim")

peaks <- readRDS(ppath)

fragments <- CreateFragmentObject(fpath)

FeatureMatrix(
  fragments = fragments,
  features = peaks
)

```

## Note

SignacSlim only provides basic matrix level processing function as well as some quality control function.

If you have any questions about how to use SignacSlim, please see [Signac](https://satijalab.org/signac/index.html) and [Seurat](https://satijalab.org/seurat/).

The copyright belongs to Tim Stuart.













