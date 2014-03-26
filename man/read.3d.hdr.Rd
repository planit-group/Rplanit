\name{read.3d.hdr}
\alias{read.3d.hdr}
\title{Read 3d array (Analyze/Gate)}
\usage{
read.3d.hdr(file.name, variable = "Dose[Gy]", voxel.origin = c(0, 0, 0))
}
\arguments{
  \item{file.name}{the file name}

  \item{variable}{name of the variable to be assigned to
  the 3d data}

  \item{voxel.origin}{coordinates of the "origin" voxel,
  i.e. the voxel cooresponding to indices (1,1,1)}
}
\description{
Read file with format "*.img, *.hdr" (Analyze/Gate) and
return a \code{values} object.
}
\seealso{
Other R/W Arrays: \code{\link{read.1d}};
\code{\link{read.3d.df}}; \code{\link{read.3d.dicom}};
\code{\link{read.3d.spliced}}; \code{\link{read.3d}}
}

