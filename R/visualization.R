#' Plot signal enrichment around TSSs
#'
#' Plot the normalized TSS enrichment score at each position relative to the
#' TSS. Requires that \code{\link{TSSEnrichment}} has already been run on the
#' assay.
#'
#' @param object A Seurat object
#' @param assay Name of the assay to use. Should have the TSS enrichment
#' information for each cell
#' already computed by running \code{\link{TSSEnrichment}}
#' @param group.by Set of identities to group cells by
#' @param idents Set of identities to include in the plot
#'
#' @importFrom SeuratObject GetAssayData DefaultAssay
#' @importFrom Matrix colMeans
#' @importFrom ggplot2 ggplot aes geom_line xlab ylab theme_classic ggtitle facet_wrap
#' theme element_blank
#'
#' @return Returns a \code{\link[ggplot2]{ggplot2}} object
#' @export
#' @concept visualization
#' @concept qc
#' @examples
#' print("see https://satijalab.org/signac/reference/tssplot")
TSSPlot <- function(
        object,
        assay = NULL,
        group.by = NULL,
        idents = NULL
) {
    assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
    if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
        stop("The requested assay is not a ChromatinAssay.")
    }
    # get the normalized TSS enrichment matrix
    positionEnrichment <- GetAssayData(
        object = object,
        assay = assay,
        slot = "positionEnrichment"
    )
    if (!("TSS" %in% names(x = positionEnrichment))) {
        stop("Position enrichment matrix not present in assay")
    }
    enrichment.matrix <- positionEnrichment[["TSS"]]

    # remove motif and expected
    if (nrow(x = enrichment.matrix) == (ncol(x = object) + 2)) {
        enrichment.matrix <- enrichment.matrix[seq((nrow(x = enrichment.matrix) - 2)), ]
    }

    # average the signal per group per base
    obj.groups <- GetGroups(
        object = object,
        group.by = group.by,
        idents = idents
    )
    groupmeans <- ApplyMatrixByGroup(
        mat = enrichment.matrix,
        groups = obj.groups,
        fun = colMeans,
        normalize = FALSE
    )

    p <- ggplot(
        data = groupmeans,
        mapping = aes(x = position, y = norm.value, color = group)
    ) +
        geom_line(stat = "identity", size = 0.2) +
        facet_wrap(facets = ~group) +
        xlab("Distance from TSS (bp)") +
        ylab(label = "Mean TSS enrichment score") +
        theme_classic() +
        theme(
            legend.position = "none",
            strip.background = element_blank()
        ) +
        ggtitle("TSS enrichment")
    return(p)
}



#' Plot fragment length histogram
#'
#' Plot the frequency that fragments of different lengths are present for
#' different groups of cells.
#'
#' @param object A Seurat object
#' @param assay Which assay to use. Default is the active assay.
#' @param region Genomic range to use. Default is fist two megabases of
#' chromosome 1. Can be a GRanges object, a string, or a vector
#' of strings.
#' @param cells Which cells to plot. Default all cells
#' @param group.by Name of one or more metadata columns to group (color) the
#' cells by. Default is the current cell identities
#' @param log.scale Display Y-axis on log scale. Default is FALSE.
#' @param ... Arguments passed to other functions
#'
#' @importFrom ggplot2 ggplot geom_histogram theme_classic aes facet_wrap xlim
#' scale_y_log10 theme element_blank
#' @importFrom SeuratObject DefaultAssay
#'
#' @export
#' @concept visualization
#' @concept qc
#' @return Returns a \code{\link[ggplot2]{ggplot}} object
#' @examples
#' print("see https://satijalab.org/signac/reference/fragmenthistogram")
FragmentHistogram <- function(
        object,
        assay = NULL,
        region = "chr1-1-2000000",
        group.by = NULL,
        cells = NULL,
        log.scale = FALSE,
        ...
) {
    cells <- SetIfNull(x = cells, y = colnames(x = object))
    assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
    if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
        stop("The requested assay is not a ChromatinAssay.")
    }
    reads <- MultiGetReadsInRegion(
        object = object,
        assay = assay,
        region = region,
        cells = cells,
        verbose = FALSE,
        ...
    )
    # add group information
    if (is.null(x = group.by)) {
        groups <- Idents(object = object)
    } else {
        md <- object[[]]
        groups <- object[[group.by]]
        groups <- groups[, 1]
        names(x = groups) <- rownames(x = md)
    }
    reads$group <- groups[reads$cell]
    if (length(x = unique(x = reads$group)) == 1) {
        p <- ggplot(data = reads, aes(length)) +
            geom_histogram(bins = 200)
    } else {
        p <- ggplot(data = reads, mapping = aes(x = length, fill = group)) +
            geom_histogram(bins = 200) +
            facet_wrap(~group, scales = "free_y")
    }
    p <- p + xlim(c(0, 800)) +
        theme_classic() +
        theme(
            legend.position = "none",
            strip.background = element_blank()
        ) +
        xlab("Fragment length (bp)") +
        ylab("Count")
    if (log.scale) {
        p <- p + scale_y_log10()
    }
    return(p)
}


















