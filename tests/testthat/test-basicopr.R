context("test basic operation of SignacSlim object")



test_that("test basic operation",{

    # get counts in peaks
    fpath <- system.file("extdata", "fragments.tsv.gz", package="SignacSlim")
    ppath <- system.file("extdata", "peaks.rds", package="SignacSlim")
    peaks <- readRDS(ppath)
    fragments <- CreateFragmentObject(fpath)
    out <- FeatureMatrix(
        fragments = fragments,
        features = peaks
    )

    expect_true(inherits(out,"dgCMatrix"))



})
