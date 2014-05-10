\name{generate.contours}
\alias{generate.contours}
\title{Generate the contours}
\usage{
generate.contours(contours.df, CT = NULL, z = NULL,
  tissue = "R3PIDEd0.2v0+MKM")
}
\arguments{
  \item{contours.df}{the dataframe in which the definition
  of the volumes is stored.}

  \item{CT}{a CT object.}

  \item{z}{a vector that specifies the z of the axial
  planes.}

  \item{tissue}{the tissue name.}
}
\value{
a contours object
}
\description{
It generates a contours object from the definition of the
VOIs. The VOIs are box-shaped an are defined in a
data.frame with columns: x_min, x_max, y_min, y_max, z_min,
z_max, voi. Each row correspond to a VOI. It needs, to
identify the slices (axial slices), a reference CT.
Alternatively the z vector can be used directly. Optionally
a tissue name can be supplied.
}
\seealso{
Other Phantoms: \code{\link{generate.ct.contours}};
\code{\link{generate.ct}}
}

