\name{display.slice.all}
\alias{display.slice.all}
\title{Display a 2D slice superimposed to the CT}
\usage{
display.slice.all(ct = NULL, contours = NULL, values = NULL,
  variable = NULL, plan = NULL, x = NA, y = NA, z = NA,
  cont = FALSE, isoc = FALSE, xlim = NULL, ylim = NULL, zlim = NULL,
  vlim = NULL, file.name = NULL, width = 7, height = 7,
  linear.graymap = TRUE, window.center = 50, window.width = 800,
  alpha.threshold = 0.7, invert.y.axis = FALSE, contour.color = "green")
}
\arguments{
  \item{ct}{the ct object (optional, if plan is specified)}

  \item{contours}{the contours object (optional)}

  \item{values}{the values object (optional, if plan is
  specified)}

  \item{variable}{the specific variable to be displayed
  (not needed if values contains only one variable, or if
  the variable is the first of the list)}

  \item{Nx,Ny,Nz}{slice index to be displayed. Only one of
  them should be defined, corresponding to an axial,
  coronal or sagittal slice (optional)}

  \item{x,y,z}{slice coordinate to be displayed. Only one
  of them should be defined, corresponding to an axial,
  coronal or sagittal slice (optional)}

  \item{plan}{the plan object (optional, if \code{values}
  and \code{ct} are defined)}

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

  \item{linear.graymap}{Uses a linear gray scale for the
  visualization of the CT Hounsfield number. If it is set
  to FALSE it uses a non linear mapping.}

  \item{window.center}{Center of the Hounsfiled window
  values to visualize. It is used if
  \code{linear.graymap=FALSE}.}

  \item{window.width}{Width of the Hounsfiled window values
  to visualize. It is used if \code{linear.graymap=FALSE}.}

  \item{alpha.threshold}{fraction of the values below which
  apply transparency. The opacity (alpha) is evaluated as a
  linear ramp (0 to 1) from 0 to \code{alpha.threshold} of
  the values range, and then alpha=1 from
  \code{alpha.threshold} to 1 of the values range.}

  \item{invert.y.axis}{Invert the y axis.}

  \item{contour.color}{The color used for the contour of
  the VOI.}
}
\description{
Display an arbitrary 2D slice (axial, coronal, sagittal) of
the values stored in a \code{values} object. The slice is
superimposed to the associated CT slice. Optionally ROI
contours can be also displayed.
}
\details{
If slice index and slice coordinate are not specified, but
the plan is specified, the displayed slice is axial and at
the isocenter of the plan. If slice index, slice coordinate
and plan are not specified, the displayed slice is axial
and at the middle of the z range.

If the plan is specified, it is possible to not specify
some or all of the following objects: values, CT and
contours. In this case the missing ones are replaced with
those associated to the plan.
}
\seealso{
Other display slices: \code{\link{display.slice.ct}};
\code{\link{display.slice}};
\code{\link{dumpslices.3d.2.png}}
}

