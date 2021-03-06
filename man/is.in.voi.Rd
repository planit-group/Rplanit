% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vois-functions.R
\name{is.in.voi}
\alias{is.in.voi}
\title{Check if point(s) is(are) inside a VOI}
\usage{
is.in.voi(x, y, z, contours, voi, progressbar = FALSE)
}
\arguments{
\item{x, y, z}{coordinate(s) of the point(s).}

\item{contours}{Contours object.}

\item{voi}{The VOI name.}
}
\value{
A logical vector.
}
\description{
Check if point(s) is(are) inside a VOI
}
\seealso{
Other VOIs: \code{\link{add.vois}},
  \code{\link{create.vois.from.logical.array}},
  \code{\link{create.vois}}, \code{\link{get.subset.vois}},
  \code{\link{get.voi.logical}},
  \code{\link{get.voiindex}}, \code{\link{get.vois}}
}
\concept{VOIs}
