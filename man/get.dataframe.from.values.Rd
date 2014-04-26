\name{get.dataframe.from.values}
\alias{get.dataframe.from.values}
\title{genera struttura dataframe da oggetto values}
\usage{
get.dataframe.from.values(values = values, vois = vois, variables = NULL,
  rois = NULL)
}
\description{
l'estrazione viene fatta sui VOI specificati ogni riga
corrisponde ad un voxel. nel dataframe ogni voxel è
associato a un VOI. voxel associati a più voi vengono
ripetuti su più righe.
}
\seealso{
Other Values: \code{\link{create.values}};
\code{\link{get.sparse.array.from.values}};
\code{\link{get.sparse.array.gate}};
\code{\link{get.values.dataframe}};
\code{\link{get.values.from.sparse.array}};
\code{\link{get.values.gate.plan}};
\code{\link{get.values.plankit.plan}};
\code{\link{get.values}}; \code{\link{merge.values}};
\code{\link{sanitize.values}}
}

