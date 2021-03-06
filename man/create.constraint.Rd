% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prescription-functions.R
\name{create.constraint}
\alias{create.constraint}
\title{Crea oggetto prescrizione (dataframe)}
\usage{
create.constraint(N = 1, VOI = "PTV", VOIIndex = NA,
  type = "EXACT", variable = "Dose[Gy]", value = 2,
  volumeFraction = 1, weight = 1)
}
\description{
nota: si tratta di un dataframe: diversi fields sono semplicemente aggiunti
appendendo delle righe al dataframe...
}
\seealso{
Other Prescription: \code{\link{check.prescription}},
  \code{\link{get.prescription}}
}
\concept{Prescription}
