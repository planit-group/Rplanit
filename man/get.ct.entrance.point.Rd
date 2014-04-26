\name{get.ct.entrance.point}
\alias{get.ct.entrance.point}
\title{Get entrance point in CT}
\usage{
get.ct.entrance.point(rays, ct)
}
\arguments{
  \item{rays}{the rays data.frame}
}
\value{
data frame of rays in wich the coordinates of the entrance
point are added (Xin, Yin, Zin).
}
\description{
Get the coordinates of the entrance point in the CT for the
specified rays.
}
\seealso{
Other Beams: \code{\link{add.field}};
\code{\link{convert.fluence2mu}};
\code{\link{convert.mu2fluence}};
\code{\link{create.beam}}; \code{\link{create.field}};
\code{\link{display.beamports}};
\code{\link{display.beams}}; \code{\link{get.beams}};
\code{\link{get.isocenter}}; \code{\link{get.rays}};
\code{\link{read.beams.fluka}}; \code{\link{read.beams}};
\code{\link{write.beams}}
}
