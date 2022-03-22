#' Get an Assay object from a given Seurat object.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns an Assay object
#'
#' @rdname GetAssay
#' @export GetAssay
#' @examples
#' if(FLASE){
#'     print("see https://satijalab.org/signac/articles/data_structures.html")
#' }
GetAssay <- function(object, ...) {
    UseMethod(generic = 'GetAssay', object = object)
}



#' Quantify aggregated genome tiles
#'
#' Quantifies fragment counts per cell in fixed-size genome bins across the
#' whole genome, then removes bins with less than a desired minimum number of
#' counts in the bin, then merges adjacent tiles into a single region.
#'
#' @param object A Seurat object or ChromatinAssay object
#' @param ... Additional arguments passed to other methods
#' @export AggregateTiles
#' @rdname AggregateTiles
#' @return AggregateTiles
#' @examples
#' if(FLASE){
#'     print("see https://satijalab.org/signac/reference/aggregatetiles")
#' }
AggregateTiles <- function(object, ...) {
    UseMethod(generic = "AggregateTiles", object = object)
}

#' Convert objects to a ChromatinAssay
#' @param x An object to convert to class \code{\link{ChromatinAssay}}
#' @param ... Arguments passed to other methods
#' @rdname as.ChromatinAssay
#' @export as.ChromatinAssay
#' @return a ChromatinAssay
as.ChromatinAssay <- function(x, ...) {
    UseMethod(generic = "as.ChromatinAssay", object = x)
}

#' Compute allele frequencies per cell
#'
#' Collapses allele counts for each strand and normalize by the total number of
#' counts at each nucleotide position.
#'
#' @param object A Seurat object, Assay, or matrix
#' @param variants A character vector of informative variants to keep. For
#' example, \code{c("627G>A","709G>A","1045G>A","1793G>A")}.
#' @param ... Arguments passed to other methods
#'
#' @export
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object with a new assay
#' containing the allele frequencies for the informative variants.
#' @examples
#' if(FLASE){
#'     print("see https://satijalab.org/signac/reference/allelefreq")
#' }
AlleleFreq <- function(object, ...) {
    UseMethod(generic = "AlleleFreq", object = object)
}

#' Annotation
#'
#' Get the annotation from a ChromatinAssay
#'
#' @param ... Arguments passed to other methods
#' @return Returns a \code{\link[GenomicRanges]{GRanges}} object
#' if the annotation data is present, otherwise returns NULL
#' @rdname Annotation
#' @export Annotation
Annotation <- function(object, ...) {
    UseMethod(generic = "Annotation", object = object)
}

#' @param value A value to set. Can be NULL, to remove the current annotation
#' information, or a \code{\link[GenomicRanges]{GRanges}} object. If a
#' \code{GRanges} object is supplied and the genome information is stored in the
#' assay, the genome of the new annotations must match the genome of the assay.
#'
#' @rdname Annotation
#' @return Returns a \code{\link[GenomicRanges]{GRanges}} object
#' if the annotation data is present, otherwise returns NULL
#' @export Annotation<-
#'
"Annotation<-" <- function(object, ..., value) {
    UseMethod(generic = 'Annotation<-', object = object)
}

#' Binarize counts
#'
#' Set counts >1 to 1 in a count matrix
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object
#' @rdname BinarizeCounts
#' @export BinarizeCounts
#' @examples
#' if(FLASE){
#'     print("see https://satijalab.org/signac/reference/binarizecounts")
#' }
BinarizeCounts <- function(object, ...) {
    UseMethod(generic = "BinarizeCounts", object = object)
}



#' Set and get cell barcode information for a Fragment object
#'
#' @param x A Seurat object
#' @param value A character vector of cell barcodes
#' @param ... Arguments passed to other methods
#' @return Cells
#' @export Cells<-
#' @examples
#' if(FLASE){
#'     print("see https://satijalab.org/seurat/articles/essential_commands.html")
#' }
"Cells<-" <- function(x, ..., value) {
    UseMethod(generic = "Cells<-", object = x)
}






#' Get the Fragment objects
#'
#' @param ... Arguments passed to other methods
#' @return Returns a list of \code{\link{Fragment}} objects. If there are
#' no Fragment objects present, returns an empty list.
#' @rdname Fragments
#' @export Fragments
Fragments <- function(object, ...) {
    UseMethod(generic = "Fragments", object = object)
}

#' @param value A \code{\link{Fragment}} object or list of Fragment objects
#' @return Fragments
#' @rdname Fragments
#' @export Fragments<-
#'
"Fragments<-" <- function(object, ..., value) {
    UseMethod(generic = 'Fragments<-', object = object)
}





#' Get or set links information
#'
#' Get or set the genomic link information for a Seurat object or ChromatinAssay
#'
#' @param ... Arguments passed to other methods
#' @return genomic link
#' @rdname Links
#' @export Links
Links <- function(object, ...) {
    UseMethod(generic = "Links", object = object)
}

#' @param value A \code{\link[GenomicRanges]{GRanges}} object
#' @rdname Links
#' @return genomic link
#' @export Links<-
"Links<-" <- function(object, ..., value) {
    UseMethod(generic = "Links<-", object = object)
}



#' Region enrichment analysis
#'
#' Count fragments within a set of regions for different groups of
#' cells.
#'
#' @param object A Seurat or ChromatinAssay object
#' @param ... Arguments passed to other methods
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object
#' @rdname RegionMatrix
#' @export RegionMatrix
#' @examples
#' if(FLASE){
#'     print("see https://satijalab.org/signac/articles/data_structures.html")
#' }

RegionMatrix <- function(object, ...) {
    UseMethod(generic = "RegionMatrix", object = object)
}

#' Compute base composition information for genomic ranges
#'
#' Compute the GC content, region lengths, and dinucleotide base frequencies
#' for regions in the assay and add to the feature metadata.
#'
#' @param object A Seurat object, Assay object, or set of genomic ranges
#' @param ... Arguments passed to other methods
#' @return Returns a dataframe
#' @rdname RegionStats
#' @export RegionStats
#' @examples
#' if(FLASE){
#'     print("see https://satijalab.org/signac/articles/data_structures.html")
#' }
RegionStats <- function(object, ...) {
    UseMethod(generic = "RegionStats", object = object)
}












