% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vois-functions.R
\name{create.vois}
\alias{create.vois}
\title{Create a VOIS object from a contours object}
\usage{
create.vois(contours, x, y, z, vois = NULL)
}
\arguments{
\item{contours}{The contours object.}

\item{x, y, z}{The coordinates of the x,y,z axes for the VOIS 3D grid.}

\item{vois}{Vois names vector. Optionally it is possible to specificy a specific set of vois to use.}
}
\value{
a VOIS object.
}
\description{
Create a VOIS object from a contours object
}
\seealso{
Other VOIs: \code{\link{add.vois}},
  \code{\link{create.vois.from.logical.array}},
  \code{\link{get.subset.vois}},
  \code{\link{get.voi.logical}},
  \code{\link{get.voiindex}}, \code{\link{get.vois}},
  \code{\link{is.in.voi}}
}
\concept{VOIs}
