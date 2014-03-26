\name{write.beamLUT}
\alias{write.beamLUT}
\title{scrive le matrici V_{b,v} (beamLUT)}
\usage{
write.beamLUT(plan, threshold = 0, threshold.variable = "Dose[Gy]",
  variable = "Dose[Gy]", file.name)
}
\description{
nota: è usata la "threshold.variable" (normalmente la dose)
per fissare la soglia le lut sono normalizzate per singolo
primario le lut sono calcolate su tutto il volume usato in
"plan"
}
\seealso{
Other BeamLUT: \code{\link{get.values.for.beam}}
}

