\name{run.gate.forward}
\alias{run.gate.forward}
\title{Evaluate forward planning (Gate)}
\usage{
run.gate.forward(plan = plan, N = plan$TotalNumberOfPrimaries,
  save.sparse.arrays = FALSE, outmessages = FALSE)
}
\arguments{
  \item{plan}{The plan object (Gate)}

  \item{N}{The mumber of primary particles to simulate
  (events)}

  \item{save.sparse.array}{Save a sparse array for the
  values (in which the information for the individual beam
  contribution is stored).}

  \item{outmessages}{Show the optuputs from the Gate
  simulation}
}
\value{
The updated plan object (Gate)
}
\description{
Evaluate forward planning (Gate)
}
\seealso{
Other PlanGate: \code{\link{create.gate.structure}};
\code{\link{create.plan.gate}};
\code{\link{create.waterbox.gate}};
\code{\link{file.variable.gate}}; \code{\link{set.ct.gate}}
}

