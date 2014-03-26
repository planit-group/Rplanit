\name{create.field}
\alias{create.field}
\title{Crea oggetto field (dataframe)}
\usage{
create.field(targetVOI = "PTV", targetVOIIndex = NA, iecGantryAngle = 0,
  iecPatientSupportAngle = 0, interSpotSpacing = c(1.5, 1.5, 1),
  spotsExtensionOutsideTarget = 0, targetIsocenter.x = NA,
  targetIsocenter.y = NA, targetIsocenter.z = NA)
}
\description{
nota: si tratta di un dataframe: diversi fields sono
semplicemente aggiunti appendendo delle righe al
dataframe...
}
\seealso{
Other Beams: \code{\link{add.field}};
\code{\link{convert.fluence2mu}};
\code{\link{convert.mu2fluence}};
\code{\link{create.beam}}; \code{\link{display.beams}};
\code{\link{get.beams}}; \code{\link{get.isocenter}};
\code{\link{read.beams}}; \code{\link{write.beams}}
}

