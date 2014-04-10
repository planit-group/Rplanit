\name{get.profile}
\alias{get.profile}
\title{Get profile}
\usage{
get.profile(values, variable = NULL, x = NULL, y = NULL, z = NULL,
  integrate = FALSE, ray = NULL, return.voxel.index = FALSE)
}
\arguments{
  \item{values}{the \code{values object}}

  \item{variable}{the variable to be profiled. If it is not
  specified, the fist available variable will be used.}

  \item{x,y,z}{definition of the axis (only two coordinates
  should be defined). For example \code{x=4.4} and
  \code{z=0.5} specifies a profile along the y-axis passing
  through x=4.4 and z=0.5}

  \item{integrate}{perform an integration over the two
  coordinates (boolean, optional). Note: integration is not
  yet implemented around an arbitrary ray.}

  \item{ray}{The ray object, Note the ray is a data.frame
  defined by 6 components (X, Y, Z, xn, yn, zn), i.e. point
  coordinates + normalized vector direction.}

  \item{return.voxel.index}{Returns the index of the voxel
  crossed by the ray.}
}
\value{
a dataframe containing the specified profile, or the
sequence of the voxel index crossed by the ray if
return.voxel.index=TRUE.
}
\description{
Get a profile (values along a specific axis, or voxel line)
from a \code{values} object. The axis is specifed trough
the definition of two of the following three coordinates:
x, y, z, or with a single ray object for arbitrary
directions.
}

