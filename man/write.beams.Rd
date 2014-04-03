\name{write.beams}
\alias{write.beams}
\title{Write beams (RTIonPlan)}
\usage{
write.beams(beams, file.name, format = "plankit", ion = "1H")
}
\arguments{
  \item{beams}{The beams object}

  \item{file.name}{The prefix of the file name.}

  \item{format}{The format. It can be 'plankit', 'gate' or
  'fluka'.}

  \item{ion.}{The primary ion. (used by fluka).}
}
\description{
Write the beams object on file. It can use different
formats (default PlanKIT).
}
\seealso{
Other Beams: \code{\link{add.field}};
\code{\link{convert.fluence2mu}};
\code{\link{convert.mu2fluence}};
\code{\link{create.beam}}; \code{\link{create.field}};
\code{\link{display.beamports}};
\code{\link{display.beams}}; \code{\link{get.beams}};
\code{\link{get.isocenter}};
\code{\link{read.beams.fluka}}; \code{\link{read.beams}}
}

