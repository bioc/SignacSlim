% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{ChunkGRanges}
\alias{ChunkGRanges}
\title{Split a genomic ranges object into evenly sized chunks}
\usage{
ChunkGRanges(granges, nchunk)
}
\arguments{
\item{granges}{A GRanges object}

\item{nchunk}{Number of chunks to split into}
}
\value{
Returns a list of GRanges objects
}
\description{
Split a genomic ranges object into evenly sized chunks
}
