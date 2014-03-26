\name{create.plan}
\alias{create.plan}
\title{Generate Plan template (PlanKIT)}
\usage{
create.plan(plan = NA, name = "test", ctFile = "ct.3d",
  contoursFile = "c1.contours",
  hounsfieldToDensityFile = "fluka-schneider.HNToDensity.1d",
  hounsfieldToStoichiometryFile = "fluka-schneider.HNToStoichiometry.1d",
  computingGridVoxelSizes = c(2, 2, 2), computingGridCoverage = 0,
  computingGridCoverageBoundary = 0, beamLine = "protonCnao",
  computingValues = c("Dose[Gy]"), cutOffRadius = 100, fields = NULL,
  prescription = NULL, maxOptimizationIterations = 50,
  outputBeamsFile = NULL, inputBeamsFile = NULL)
}
\value{
A Plan object (list).
}
\description{
Generate a default (empty) template of the Plan to be used
with PlanKIT.
}
\seealso{
Other Plan: \code{\link{read.plan}};
\code{\link{save.plan}}
}

