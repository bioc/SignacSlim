% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{TSSPlot}
\alias{TSSPlot}
\title{Plot signal enrichment around TSSs}
\usage{
TSSPlot(object, assay = NULL, group.by = NULL, idents = NULL)
}
\arguments{
\item{object}{A Seurat object}

\item{assay}{Name of the assay to use. Should have the TSS enrichment
information for each cell
already computed by running \code{\link{TSSEnrichment}}}

\item{group.by}{Set of identities to group cells by}

\item{idents}{Set of identities to include in the plot}
}
\value{
Returns a \code{\link[ggplot2]{ggplot2}} object
}
\description{
Plot the normalized TSS enrichment score at each position relative to the
TSS. Requires that \code{\link{TSSEnrichment}} has already been run on the
assay.
}
\examples{
print("see https://satijalab.org/signac/reference/tssplot")
}
\concept{qc}
\concept{visualization}
