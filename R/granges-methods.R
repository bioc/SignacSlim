#' @include generics.R
#' @importFrom SeuratObject DefaultAssay
#' @importFrom GenomicRanges granges
NULL

#' Access genomic ranges for ChromatinAssay objects
#'
#' Methods for accessing \code{\link[GenomicRanges]{GRanges}} object
#' information stored in a \code{\link{ChromatinAssay}} object.
#'
#' @name granges-methods
#' @param x A \code{\link{ChromatinAssay}} object
#' @param use.names Whether the names on the genomic ranges should be
#' propagated to the returned object.
#' @param use.mcols Not supported for \code{\link{ChromatinAssay}} objects
#' @param ... Additional arguments
#'
#' @return Returns a \code{\link[GenomicRanges]{GRanges}} object
#'
#' @aliases granges granges,ChromatinAssay-method
#' @seealso
#' \itemize{
#'   \item{\link[GenomicRanges]{granges} in the \pkg{GenomicRanges} package.}
#'   \item{\link{ChromatinAssay-class}}
#'  }
#' @importFrom GenomicRanges granges
#' @exportMethod granges
#' @concept granges
#' @examples
#' ppath <- system.file("extdata", "peaks.rds", package="SignacSlim")
#' peaks <- readRDS(ppath)
setMethod(
  f = "granges",
  signature = "ChromatinAssay",
  definition = function(x, use.names = TRUE, use.mcols = FALSE, ...) {
    if (!identical(x = use.mcols, y = FALSE)) {
      stop("\"granges\" method for ChromatinAssay objects ",
           "does not support the 'use.mcols' argument")
    }
    slot(object = x, name = "ranges")
  }
)

