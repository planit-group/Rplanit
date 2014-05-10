\name{get.values.gate.plan}
\alias{get.values.gate.plan}
\title{Get values object (Gate)}
\usage{
\method{get.values}{gate.plan}(plan, center.coordinates = NULL)
}
\arguments{
  \item{plan.gate}{The plan object.}

  \item{center.coordinates}{The coordinates of the center
  of the values 3D distribution. If it is NULL, the center
  of the CT is assumed.}
}
\description{
Get the values matrices for a plan (Gate).
}
\seealso{
Other Values: \code{\link{create.values}};
\code{\link{get.dataframe.from.values}};
\code{\link{get.sparse.array.from.values}};
\code{\link{get.sparse.array.gate}};
\code{\link{get.values.dataframe}};
\code{\link{get.values.from.sparse.array}};
\code{\link{get.values.plankit.plan}};
\code{\link{get.values}}; \code{\link{merge.values}};
\code{\link{sanitize.values}}
}

