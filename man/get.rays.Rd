\name{get.rays}
\alias{get.rays}
\title{Get rays from beams}
\usage{
get.rays(beams, s = 4000, unique = FALSE)
}
\arguments{
  \item{beams}{the beams data.frame}

  \item{s}{distance of the source from the isocenter (it
  uses only one virtual source).}
}
\description{
Get rays (normalized vector directions + absolute
coordinate) in the reference frame of the CT.
}
\seealso{
Other Beams: \code{\link{add.field}};
\code{\link{convert.fluence2mu}};
\code{\link{convert.mu2fluence}};
\code{\link{create.beam}}; \code{\link{create.field}};
\code{\link{display.beamports}};
\code{\link{display.beams}}; \code{\link{get.beams}};
\code{\link{get.isocenter}};
\code{\link{read.beams.fluka}}; \code{\link{read.beams}};
\code{\link{write.beams}}
}

