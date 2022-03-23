#' Construct a feature x cell matrix from a genomic fragments file
#'
#' @param fragments A list of \code{\link{Fragment}} objects. Note that if
#' setting the \code{cells} parameter, the requested cells should be present in
#' the supplied \code{Fragment} objects. However, if the cells information in
#' the fragment object is not set (\code{Cells(fragments)} is \code{NULL}), then
#' the fragment object will still be searched.
#' @param features A GRanges object containing a set of genomic intervals.
#' These will form the rows of the matrix, with each entry recording the number
#' of unique reads falling in the genomic region for each cell.
#' @param cells Vector of cells to include. If NULL, include all cells found
#' in the fragments file
#' @param process_n Number of regions to load into memory at a time, per thread.
#' Processing more regions at once can be faster but uses more memory.
#' @param sep Vector of separators to use for genomic string. First element is
#' used to separate chromosome and coordinates, second separator is used to
#' separate start and end coordinates.
#' @param verbose Display messages
#'
#' @return Returns a sparse matrix
#'
#' @export
#' @importFrom SeuratObject RowMergeSparseMatrices
#'
#' @examples
#' \dontrun{
#' # test data from package Signac
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' fragments <- CreateFragmentObject(fpath)
#' FeatureMatrix(
#'   fragments = fragments,
#'   features = granges(atac_small)
#' )
#' }
#'
FeatureMatrix <- function(
        fragments,
        features,
        cells = NULL,
        process_n = 2000,
        sep = c("-", "-"),
        verbose = TRUE
) {
    if (!inherits(x = features, what = "GRanges")) {
        if (inherits(x = features, what = "character")) {
            features <- StringToGRanges(
                regions = features,
                sep = sep
            )
        } else {
            stop("features should be a GRanges object")
        }
    }
    if (!inherits(x = fragments, what = "list")) {
        if (inherits(x = fragments, what = "Fragment")) {
            fragments <- list(fragments)
        } else {
            stop("fragments should be a list of Fragment objects")
        }
    }
    # if cells is not NULL, iterate over all fragment objects
    # and find which objects contain cells that are requested
    if (!is.null(x = cells)) {
        obj.use <- c()
        for (i in seq_along(along.with = fragments)) {
            if (is.null(x = Cells(fragments[[i]]))) {
                # cells information not set for fragment object
                obj.use <- c(obj.use, i)
            } else if (any(cells %in% Cells(x = fragments[[i]]))) {
                obj.use <- c(obj.use, i)
            }
        }
    } else {
        obj.use <- seq_along(along.with = fragments)
    }
    # create a matrix from each fragment file
    mat.list <- sapply(
        X = obj.use,
        FUN = function(x) {
            SingleFeatureMatrix(
                fragment = fragments[[x]],
                features = features,
                cells = cells,
                sep = sep,
                verbose = verbose,
                process_n = process_n
            )
        })
    # merge all the matrices
    if (length(x = mat.list) == 1) {
        return(mat.list[[1]])
    } else {
        featmat <- Reduce(f = RowMergeSparseMatrices, x = mat.list)
        return(featmat)
    }
}
