% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vois-functions.R
\name{get.contours}
\alias{get.contours}
\title{Get contours from plan}
\usage{
get.contours(plan)
}
\arguments{
\item{plan}{the \code{plan} object. It can be a list of plans.}
}
\value{
A \code{contours} dataframe (or a list of dataframes) consisting of:
\item{id}{index of the contours}
\item{polygon}{index of the polygon}
\item{slice}{index of the slice}
\item{x,y,z}{x,y,z coordinates}
\item{contour}{contour (VOI) name}
\item{tissue}{tissue+biological model name}
\item{type}{contour type (target, OAR, etc.)}
}
\description{
Get the contours dataframe from the \code{plan} object.
}
\seealso{
Other Contours: \code{\link{add.colours.contours}},
  \code{\link{display.contours.legend}},
  \code{\link{sanitize.contours}},
  \code{\link{volume.contours}}
}
\concept{Contours}
