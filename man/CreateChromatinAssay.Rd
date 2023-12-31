% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\name{CreateChromatinAssay}
\alias{CreateChromatinAssay}
\title{Create ChromatinAssay object}
\usage{
CreateChromatinAssay(
  counts,
  data,
  min.cells = 0,
  min.features = 0,
  max.cells = NULL,
  ranges = NULL,
  motifs = NULL,
  fragments = NULL,
  genome = NULL,
  annotation = NULL,
  bias = NULL,
  positionEnrichment = NULL,
  sep = c("-", "-"),
  validate.fragments = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{counts}{Unnormalized data (raw counts)}

\item{data}{Normalized data; if provided, do not pass counts}

\item{min.cells}{Include features detected in at least this many cells.
Will subset the counts matrix as well.
To reintroduce excluded features, create a new object with a lower cutoff.}

\item{min.features}{Include cells where at least this many features are
detected.}

\item{max.cells}{Include features detected in less than this many cells.
Will subset the counts matrix as well.
To reintroduce excluded features, create a new object with a higher cutoff.
This can be useful for chromatin assays where certain artefactual loci
accumulate reads in all cells. A percentage cutoff can also be set using
'q' followed by the percentage of cells, for example 'q90' will discard
features detected in 90 percent of cells.
If NULL (default), do not apply any maximum value.}

\item{ranges}{A set of \code{\link[GenomicRanges]{GRanges}} corresponding to
the rows of the input matrix}

\item{motifs}{A Motif object (not required)}

\item{fragments}{Path to a tabix-indexed fragments file for the data
contained in the input matrix. If multiple fragment files are required,
you can add additional \code{\link{Fragment}} object to the assay after it is
created using the \code{\link{CreateFragmentObject}} and
\code{\link{Fragments}} functions. Alternatively, a list of
\code{\link{Fragment}} objects can be provided.}

\item{genome}{A \code{\link[GenomeInfoDb]{Seqinfo}} object containing basic
information about the genome used. Alternatively, the name of a UCSC genome
can be provided and the sequence information will be downloaded from UCSC.}

\item{annotation}{A set of \code{\link[GenomicRanges]{GRanges}} containing
annotations for the genome used}

\item{bias}{A Tn5 integration bias matrix}

\item{positionEnrichment}{A named list of matrices containing positional
signal enrichment information for each cell. Should be a cell x position
matrix, centered on an element of interest (for example, TSS sites).}

\item{sep}{Separators to use for strings encoding genomic coordinates.
First element is used to separate the chromosome from the coordinates,
second element is used to separate the start from end coordinate. Only
used if \code{ranges} is NULL.}

\item{validate.fragments}{Check that cells in the assay are present in the
fragment file.}

\item{verbose}{Display messages}

\item{...}{Additional arguments passed to \code{\link{CreateFragmentObject}}}
}
\value{
ChromatinAssay
}
\description{
Create a \code{\link{ChromatinAssay}} object from a count matrix or
normalized data matrix. The expected format of the input matrix is features x
cells. A set of genomic ranges must be supplied along with the matrix, with
the length of the ranges equal to the number of rows in the matrix. If a set
of genomic ranges are not supplied, they will be extracted from the
row names of the matrix.
}
\examples{
print("see https://satijalab.org/signac/reference/createchromatinassay")
}
\concept{assay}
