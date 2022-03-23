#' NucleosomeSignal
#'
#' Calculate the strength of the nucleosome signal per cell.
#' Computes the ratio of fragments between 147 bp and 294 bp (mononucleosome) to
#' fragments < 147 bp (nucleosome-free)
#'
#' @param object A Seurat object
#' @param assay Name of assay to use. Only required if a fragment path is not
#' provided. If NULL, use the active assay.
#' @param n Number of lines to read from the fragment file. If NULL, read all
#' lines. Default scales with the number of cells in the object.
#' @param verbose Display messages
#' @param ... Arguments passed to other functions
#'
#' @importFrom dplyr group_by summarize
#' @importFrom stats ecdf
#'
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object with
#' added metadata for the ratio of mononucleosomal to nucleosome-free fragments
#' per cell, and the percentile rank of each ratio.
#' @export
#' @concept qc
#' @importFrom fastmatch fmatch
#' @importFrom SeuratObject AddMetaData
#' @importFrom stats ecdf
#'
#' @examples
#' print("see https://satijalab.org/signac/reference/nucleosomesignal")
NucleosomeSignal <- function(
        object,
        assay = NULL,
        n = ncol(object) * 5e3,
        verbose = TRUE,
        ...
) {
    assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
    if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
        stop("The requested assay is not a ChromatinAssay")
    }
    # first check that fragments are present
    frags <- Fragments(object = object[[assay]])
    if (length(x = frags) == 0) {
        stop("No fragment files present in assay")
    }
    verbose <- as.logical(x = verbose)
    af <- list()
    for (i in seq_along(along.with = frags)) {
        counts <- ExtractFragments(
            fragments = frags[[i]],
            n = n,
            verbose = verbose
        )
        cells.keep <- fmatch(
            x = counts$CB, table = colnames(x = object), nomatch = 0L
        )
        rownames(x = counts) <- counts$CB
        counts <- counts[
            cells.keep > 0, c("mononucleosomal", "nucleosome_free")
        ]
        af[[i]] <- counts
    }
    af <- do.call(what = rbind, args = af)
    af$nucleosome_signal <- af$mononucleosomal / af$nucleosome_free
    e.dist <- ecdf(x = af$nucleosome_signal)
    af$nucleosome_percentile <- round(
        x = e.dist(af$nucleosome_signal),
        digits = 2
    )
    af <- af[, c("nucleosome_signal", "nucleosome_percentile")]
    object <- AddMetaData(object = object, metadata = af)
    return(object)
}






























