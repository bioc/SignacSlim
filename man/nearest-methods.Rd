% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iranges-methods.R
\name{nearest-methods}
\alias{nearest-methods}
\alias{precede}
\alias{precede,ANY,ChromatinAssay-method}
\alias{precede,ChromatinAssay,ANY-method}
\alias{precede,ChromatinAssay,ChromatinAssay-method}
\alias{follow,ANY,ChromatinAssay-method}
\alias{follow}
\alias{follow,ChromatinAssay,ANY-method}
\alias{follow,ChromatinAssay,ChromatinAssay-method}
\alias{nearest,ANY,ChromatinAssay-method}
\alias{nearest}
\alias{nearest,ChromatinAssay,ANY-method}
\alias{nearest,ChromatinAssay,ChromatinAssay-method}
\alias{distance,ANY,ChromatinAssay-method}
\alias{distance}
\alias{distance,ChromatinAssay,ANY-method}
\alias{distance,ChromatinAssay,ChromatinAssay-method}
\alias{distanceToNearest,ANY,ChromatinAssay-method}
\alias{distanceToNearest}
\alias{distanceToNearest,ChromatinAssay,ANY-method}
\alias{distanceToNearest,ChromatinAssay,ChromatinAssay-method}
\title{Find the nearest range neighbors for ChromatinAssay objects}
\usage{
\S4method{precede}{ANY,ChromatinAssay}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{precede}{ChromatinAssay,ANY}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{precede}{ChromatinAssay,ChromatinAssay}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{follow}{ANY,ChromatinAssay}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{follow}{ChromatinAssay,ANY}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{follow}{ChromatinAssay,ChromatinAssay}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{nearest}{ANY,ChromatinAssay}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{nearest}{ChromatinAssay,ANY}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{nearest}{ChromatinAssay,ChromatinAssay}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{distance}{ANY,ChromatinAssay}(x, y, ignore.strand = FALSE, ...)

\S4method{distance}{ChromatinAssay,ANY}(x, y, ignore.strand = FALSE, ...)

\S4method{distance}{ChromatinAssay,ChromatinAssay}(x, y, ignore.strand = FALSE, ...)

\S4method{distanceToNearest}{ANY,ChromatinAssay}(x, subject, ignore.strand = FALSE, ...)

\S4method{distanceToNearest}{ChromatinAssay,ANY}(x, subject, ignore.strand = FALSE, ...)

\S4method{distanceToNearest}{ChromatinAssay,ChromatinAssay}(x, subject, ignore.strand = FALSE, ...)
}
\arguments{
\item{x}{A query \code{\link{ChromatinAssay}} object}

\item{subject}{The subject \code{\link[GenomicRanges]{GRanges}} or
\code{\link{ChromatinAssay}} object. If missing, \code{x} is used as the
subject.}

\item{select}{Logic for handling ties.
See \code{\link[GenomicRanges]{nearest-methods}} in the \pkg{GenomicRanges}
package.}

\item{ignore.strand}{Logical argument controlling whether strand information
should be ignored.}

\item{y}{For the \code{distance} method, a
\code{\link[GenomicRanges]{GRanges}} object or a \code{\link{ChromatinAssay}}
object}

\item{...}{Additional arguments for methods}
}
\value{
nearest range neighbors
}
\description{
The \code{precede, follow, nearest, distance, distanceToNearest} methods
are available for \code{\link{ChromatinAssay}} objects.
}
\section{Functions}{
\itemize{
\item \code{precede,ChromatinAssay,ANY-method}: method for ChromatinAssay, ANY

\item \code{precede,ChromatinAssay,ChromatinAssay-method}: method for ChromatinAssay, ChromatinAssay

\item \code{follow,ANY,ChromatinAssay-method}: method for ANY, ChromatinAssay

\item \code{follow,ChromatinAssay,ANY-method}: method for ChromatinAssay, ANY

\item \code{follow,ChromatinAssay,ChromatinAssay-method}: method for ChromatinAssay, ChromatinAssay

\item \code{nearest,ANY,ChromatinAssay-method}: method for ANY, ChromatinAssay

\item \code{nearest,ChromatinAssay,ANY-method}: method for ChromatinAssay, ANY

\item \code{nearest,ChromatinAssay,ChromatinAssay-method}: method for ChromatinAssay, ChromatinAssay

\item \code{distance,ANY,ChromatinAssay-method}: method for ANY, ChromatinAssay

\item \code{distance,ChromatinAssay,ANY-method}: method for ChromatinAssay, ANY

\item \code{distance,ChromatinAssay,ChromatinAssay-method}: method for ChromatinAssay, ChromatinAssay

\item \code{distanceToNearest,ANY,ChromatinAssay-method}: method for ANY, ChromatinAssay

\item \code{distanceToNearest,ChromatinAssay,ANY-method}: method for ChromatinAssay, ANY

\item \code{distanceToNearest,ChromatinAssay,ChromatinAssay-method}: method for ChromatinAssay, ChromatinAssay
}}

\examples{
print("see https://satijalab.org/signac/reference/chromatinassay-class")
}
\seealso{
\itemize{
  \item{\link[IRanges]{nearest-methods} in the \pkg{IRanges} package.}
  \item{\link[GenomicRanges]{nearest-methods} in the \pkg{GenomicRanges}
  package}
  \item{\link{ChromatinAssay-class}}
 }
}
\concept{nearest}
