% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NucleosomeSignal.R
\name{NucleosomeSignal}
\alias{NucleosomeSignal}
\title{NucleosomeSignal}
\usage{
NucleosomeSignal(
  object,
  assay = NULL,
  n = ncol(object) * 5000,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{object}{A Seurat object}

\item{assay}{Name of assay to use. Only required if a fragment path is not
provided. If NULL, use the active assay.}

\item{n}{Number of lines to read from the fragment file. If NULL, read all
lines. Default scales with the number of cells in the object.}

\item{verbose}{Display messages}

\item{...}{Arguments passed to other functions}
}
\value{
Returns a \code{\link[SeuratObject]{Seurat}} object with
added metadata for the ratio of mononucleosomal to nucleosome-free fragments
per cell, and the percentile rank of each ratio.
}
\description{
Calculate the strength of the nucleosome signal per cell.
Computes the ratio of fragments between 147 bp and 294 bp (mononucleosome) to
fragments < 147 bp (nucleosome-free)
}
\examples{
print("see https://satijalab.org/signac/reference/nucleosomesignal")
}
\concept{qc}
