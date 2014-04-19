\name{create.values}
\alias{create.values}
\title{Create values object}
\usage{
create.values(array.values = NULL, variables = "variable name", x, y, z)
}
\arguments{
  \item{array}{the 3D array}

  \item{variables}{The name of the variable (can be a
  vector of names).}

  \item{x,y,z}{coordinate vectors}
}
\value{
a \code{values} object
}
\description{
Create values object from array. Note: the value
distributions are stored in a multidimensional array:
values$values with dim(values$values) = c(Nv, Nx, Nz, Ny).
}
\seealso{
Other Values: \code{\link{get.dataframe.from.values}};
\code{\link{get.sparse.array.from.values}};
\code{\link{get.sparse.array.gate}};
\code{\link{get.values.dataframe}};
\code{\link{get.values.from.sparse.array}};
\code{\link{get.values.gate}}; \code{\link{get.values}};
\code{\link{merge.values}}; \code{\link{sanitize.values}}
}

