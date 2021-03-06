% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lut-functions.R
\name{get.dose.beamLUT}
\alias{get.dose.beamLUT}
\title{Evaluate net Dose from beamLUTs  (pure-dek)}
\usage{
get.dose.beamLUT(beamLUTs, beams, x = NULL, y = NULL, z = NULL,
  variable = "Dose[Gy]")
}
\arguments{
\item{beamLUTs}{the beamLUTs data frame.}

\item{beams}{the beams data frame.}

\item{x, y, z}{the coordinates of the 3D array.}

\item{variable}{to override the default name of the variable ('Dose[Gy]').}
}
\value{
a values object (if  x, y and z are specified), otherwise a sparse array.
}
\description{
Get the net Dose distribution (values object or sparse array) from the beamLUTs. If coordinates x, y and z are specified it returns a complete values object
}
\seealso{
Other BeamLUT: \code{\link{combine.beamLUTs}},
  \code{\link{get.alpha.beta.beamLUT}},
  \code{\link{get.bdose.beamLUT}},
  \code{\link{get.beamLUTs}},
  \code{\link{get.letd.beamLUT}},
  \code{\link{get.rbe.beamLUT}},
  \code{\link{get.values.for.beam}},
  \code{\link{read.beamLUTs}}, \code{\link{read.beamLUT}},
  \code{\link{write.beamLUT}}
}
\concept{BeamLUT}
