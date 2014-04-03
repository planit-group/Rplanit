\name{display.beamports}
\alias{display.beamports}
\title{Beam-port splot}
\usage{
display.beamports(beams, plan = NULL, show.plot = TRUE, file.name = NULL,
  height = 7, width = 7)
}
\arguments{
  \item{beams}{beams dataframe.}

  \item{plan}{The plan object.}

  \item{shot.plot}{Show plot.}

  \item{file.name}{File name for saving the plot.}

  \item{height}{The height of the saved plot image
  (inches).}

  \item{width}{The width of the saved plot image (inches).}
}
\value{
If show.plot=FALSE, it returns a ggplot object.
}
\description{
Beam-port splot
}
\seealso{
Other Beams: \code{\link{add.field}};
\code{\link{convert.fluence2mu}};
\code{\link{convert.mu2fluence}};
\code{\link{create.beam}}; \code{\link{create.field}};
\code{\link{display.beams}}; \code{\link{get.beams}};
\code{\link{get.isocenter}};
\code{\link{read.beams.fluka}}; \code{\link{read.beams}};
\code{\link{write.beams}}
}

