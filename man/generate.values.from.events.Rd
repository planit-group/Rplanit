\name{generate.values.from.events}
\alias{generate.values.from.events}
\title{Create values form events}
\usage{
generate.values.from.events(Xe, Ye, Ze, xbin, ybin, zbin, weight = 1,
  variable = "Activity", sparse.matrix = FALSE)
}
\arguments{
  \item{Xe}{event x coorditate vector.}

  \item{Ye}{event y coordinate vector.}

  \item{Ze}{event z coordinate vector.}

  \item{weight}{weight vector.}

  \item{xbin}{bins x coordinates (extremes).}

  \item{ybin}{bins y coordinates (extremes).}

  \item{zbin}{bins z coordinates (extremes).}

  \item{variable}{The variable name.}

  \item{sparse.matrix}{Return a sparse matrix if TRUE.}
}
\value{
Values object or a sparse matrix.
}
\description{
Uses the (x,y,z) coordinates of a shower of events to
create a values object with the 3D histogram of the events.
}
\seealso{
Other Values, Values Utilities:
\code{\link{create.intervals}};
\code{\link{get.coordinates}};
\code{\link{get.voxel.index}}
}

