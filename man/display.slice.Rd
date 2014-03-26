\name{display.slice}
\alias{display.slice}
\title{Display a 2D slice}
\usage{
display.slice(values = NULL, variable = NULL, Nx = NA, Ny = NA,
  Nz = NA, x = NA, y = NA, z = NA, plan = NULL, bw = FALSE,
  cont = FALSE, isoc = FALSE, xlim = NULL, ylim = NULL, zlim = NULL,
  vlim = NULL, file.name = NULL, width = 7, height = 7)
}
\arguments{
  \item{values}{the values object (optional if plan is
  defined)}

  \item{variable}{the specific variable to be displayed
  (not needed if values contains only one variable, or if
  the variable is the first of the list)}

  \item{Nx,Ny,Nz}{slice index to be displayed. Only one of
  them should be defined, corresponding to an axial,
  coronal or sagittal slice (optional)}

  \item{x,y,z}{slice coordinate to be displayed. Only one
  of them should be defined, corresponding to an axial,
  coronal or sagittal slice (optional)}

  \item{plan}{the plan object (optional if values is
  defined)}

  \item{cont}{display isolevel contours of the values
  (optional, boolean)}

  \item{isoc}{display the isocenter on the slice (optional,
  boolean)}

  \item{bw}{display using gray color levels. It colud be
  used to display CT values (optional, boolean)}

  \item{xlim,ylim,zlim}{two element vectors containing the
  coordinate ranges to be displayed (optional)}

  \item{vlim}{two element vector containing the value range
  to be displayed. Values outside the range are displayed
  as 'blank', i.e. no color tile (optional)}

  \item{file.name}{name of the file on which to save the
  figure. By default no file is written (optional)}

  \item{width,height}{sizes of the figure to be saved
  (inches) (optional)}
}
\value{
used slice (index and coordinates) data.frame. It is used
to make further combination with other plots and procedures
(e.g. \code{display.slice.all})
}
\description{
Display an arbitrary 2D slice (axial, coronal, sagittal) of
the values stored in a \code{values} object.
}
\details{
If slice index and slice coordinate are not specified, but
the plan is specified, the displayed slice is axial and at
the isocenter of the plan. If slice index, slice coordinate
and plan are not specified, the displayed slice is axial
and at the middle of the z range.

If the plan is specified, it is possible to not specify the
values object. In this case the one associated to the plan
is used.
}
\seealso{
Other display slices: \code{\link{display.slice.all}};
\code{\link{display.slice.ct}};
\code{\link{dumpslices.3d.2.png}}
}

