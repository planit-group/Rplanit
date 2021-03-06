% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/display-functions.R
\name{display.contours.legend}
\alias{display.contours.legend}
\title{Add a legend with contours names/colours.}
\usage{
display.contours.legend(contours, add = TRUE, position = "topleft",
  cex = 0.8, ...)
}
\arguments{
\item{contours}{The contours dataframe.}

\item{add}{Put the legend in the active plot. If FALSE, it returns a stand alone legend in a new display.}

\item{position}{The position of the legend ('topleft', 'bottomleft', etc.)}

\item{cex}{dimension of the text font}

\item{...}{pther parameters to legend()}
}
\description{
If colours are not specitied in the contours dataframe, they will be added using the "raimbow" colour selection.
}
\seealso{
Other Contours: \code{\link{add.colours.contours}},
  \code{\link{get.contours}},
  \code{\link{sanitize.contours}},
  \code{\link{volume.contours}}
}
\concept{Contours}
