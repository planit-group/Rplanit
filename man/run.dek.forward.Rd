\name{run.dek.forward}
\alias{run.dek.forward}
\title{Forward planning (PlanKIT)}
\usage{
run.dek.forward(plan, outmessages = TRUE)
}
\arguments{
  \item{plan}{the Plan object. It should include also the
  specification of the beams (i.e. the "RTIonPlan").}

  \item{outmessages}{if TRUE, it displays the messages from
  the PlanKIT output.}
}
\value{
a Plan Object (including data/information of the evaluated
output of the forward planning, e.g. values, vois, etc.)
}
\description{
Perform a forward planning with PlanKIT.
}
\seealso{
Other PlanKIT: \code{\link{run.dek.inverse}}
}

