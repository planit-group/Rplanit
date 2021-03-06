% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reading-writing-functions.R
\name{read.3d}
\alias{read.3d}
\title{Read 3d array (PlanKIT)}
\usage{
read.3d(file.name, x.min = -Inf, x.max = +Inf, y.min = -Inf,
  y.max = +Inf, z.min = -Inf, z.max = +Inf)
}
\arguments{
\item{file.name}{the file name}

\item{x.min, x.max, y.min, y.max, z.min, z.max}{load only the volume inside the specified ranges (optional)}
}
\description{
Read file with format "3d" (PlanKIT) and return a \code{values} object.
}
\seealso{
Other R/W Arrays: \code{\link{read.1d}},
  \code{\link{read.3d.hdr}}
}
\concept{R/W Arrays}
