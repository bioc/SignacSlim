% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{Read10X_h5}
\alias{Read10X_h5}
\title{Read 10X hdf5 file}
\usage{
Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
}
\arguments{
\item{filename}{Path to h5 file}

\item{use.names}{Label row names with feature names rather than ID numbers.}

\item{unique.features}{Make feature names unique (default TRUE)}
}
\value{
Returns a sparse matrix with rows and columns labeled. If multiple
genomes are present, returns a list of sparse matrices (one per genome).
}
\description{
Read count matrix from 10X CellRanger hdf5 file.
This can be used to read both scATAC-seq and scRNA-seq matrices.
}
\examples{
print("see https://satijalab.org/seurat/reference/read10x_h5")
}
\concept{preprocessing}
