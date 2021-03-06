% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lut-functions.R
\name{get.beamLUTs}
\alias{get.beamLUTs}
\title{Get beamLUTs}
\usage{
get.beamLUTs(plan, concatenate = FALSE)
}
\arguments{
\item{plan}{The plan object.}

\item{concatenate}{use alternative faster reading approach (parameter passed to read.beamLUT(), may not work on windows and osx.)}
}
\description{
Get the evaluated beamLUTs from the plan
}
\seealso{
Other BeamLUT: \code{\link{combine.beamLUTs}},
  \code{\link{get.alpha.beta.beamLUT}},
  \code{\link{get.bdose.beamLUT}},
  \code{\link{get.dose.beamLUT}},
  \code{\link{get.letd.beamLUT}},
  \code{\link{get.rbe.beamLUT}},
  \code{\link{get.values.for.beam}},
  \code{\link{read.beamLUTs}}, \code{\link{read.beamLUT}},
  \code{\link{write.beamLUT}}
}
\concept{BeamLUT}
