\name{get.values.for.beam}
\alias{get.values.for.beam}
\title{Evaluate BeamLUT}
\usage{
get.values.for.beam(beam.index, fluence = NULL, plan,
  variables = "Dose[Gy]", temp.name = "temp", remove.temp = TRUE)
}
\arguments{
  \item{beam.index}{index of the beam (it can be a vector
  if indices)}

  \item{fluence}{the fluence (number of particles) for the
  beam(s). It can be a vector}

  \item{plan}{the plan object}

  \item{variable}{string specifing the variable to be
  evaluated. It can be a vecor of variables}

  \item{temp.name}{temporary prefix of the plan name used
  during the evaluation}

  \item{remove.temp}{TRUE to remove temporary files after
  the evaluation}
}
\value{
a values object
}
\description{
Evaluate the beamLUT of the selected beam(s) for the
specific plan. The beam(s) can be one or more (by using a
vector of indices). If the beams are more than one, the
corresponding beamLUT is given by the net contribution of
all the specified beams Optionally it is possible to
specify also the fluences (i.e. the number of particles)
for each beam.
}
\seealso{
Other BeamLUT: \code{\link{write.beamLUT}}
}

