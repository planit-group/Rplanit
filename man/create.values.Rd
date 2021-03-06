% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/values-functions.R
\name{create.values}
\alias{create.values}
\title{Create values object}
\usage{
create.values(array.values = NULL, variables = "variable name", x, y,
  z)
}
\arguments{
\item{variables}{The name of the variable (can be a vector of names).}

\item{x, y, z}{coordinate vectors}

\item{array}{the 3D array. Optional. If it is not specified an array of zeros is assumed.}
}
\value{
a \code{values} object
}
\description{
Create values object from array. Note: the value distributions are stored in a multidimensional array:
values$values with dim(values$values) = c(Nv, Nx, Nz, Ny).
}
\seealso{
Other Values: \code{\link{combine.values}},
  \code{\link{dataframe.from.values}},
  \code{\link{get.sparse.array.gate}},
  \code{\link{get.values.dataframe}},
  \code{\link{get.values.gate.plan}},
  \code{\link{get.values.list}},
  \code{\link{get.values.plankit.plan}},
  \code{\link{get.values}}, \code{\link{merge.values}},
  \code{\link{sanitize.values}},
  \code{\link{sparse.array.from.values}},
  \code{\link{threshold.values}},
  \code{\link{values.from.sparse.array}}
}
\concept{Values}
