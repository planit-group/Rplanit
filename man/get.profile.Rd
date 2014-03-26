\name{get.profile}
\alias{get.profile}
\title{Get profile}
\usage{
get.profile(values, variable = NULL, x = NULL, y = NULL, z = NULL,
  integrate = FALSE)
}
\arguments{
  \item{values}{the \code{values object}}

  \item{variable}{the variable to be profiled}

  \item{x,y,z}{definition of the axis (only two coordinates
  should be defined). For example \code{x=4.4} and
  \code{z=0.5} specifies a profile along the y-axis passing
  through x=4.4 and z=0.5}

  \item{integrate}{perform an integration over the two
  coordinates (boolean, optional)}
}
\value{
a dataframe containing the specified profile
}
\description{
Get a profile (values along a specific axis, or voxel line)
from a \code{values} object. The axis is specifed trough
the definition of two of the following three coordinates:
x, y, z.
}
\seealso{
Other Values: \code{\link{create.values}};
\code{\link{get.dataframe.from.values}};
\code{\link{get.sparse.array.from.values}};
\code{\link{get.values.dataframe}};
\code{\link{get.values.from.sparse.array}};
\code{\link{get.values}}; \code{\link{sanitize.values}}
}

