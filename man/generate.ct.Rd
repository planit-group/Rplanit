\name{generate.ct}
\alias{generate.ct}
\title{Generate a CT phantom}
\usage{
generate.ct(CT.df, deltaX, deltaY, deltaZ)
}
\arguments{
  \item{CT.df}{the dataframe in which are defined the
  volumes.}

  \item{deltaX,deltaY,deltaZ}{are the dimensions of the
  voxels.}
}
\value{
a CT (values) object
}
\description{
It uses the definition of box-shaped volumes stored in a
data.frame with columns: x_min, x_max, y_min, y_max, z_min,
z_max, HU, voi.
}
\seealso{
Other Phantoms: \code{\link{generate.contours}};
\code{\link{generate.ct.contours}}
}

