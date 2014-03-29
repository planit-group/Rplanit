\name{run.dek.inverse}
\alias{run.dek.inverse}
\title{Inverse planning (PlanKIT)}
\usage{
run.dek.inverse(plan, outmessages = FALSE)
}
\arguments{
  \item{plan}{the Plan object.}

  \item{outmessages}{if TRUE, it displays the messages from
  PlanKIT output.}
}
\value{
a Plan Object (including data/information of the evaluated
output of the inverse planning, e.g. beams, values, vois,
etc.)
}
\description{
Perform an inverse planning with PlanKIT.
}
\seealso{
Other PlanKIT: \code{\link{run.dek.forward}}
}

