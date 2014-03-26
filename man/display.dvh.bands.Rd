\name{display.dvh.bands}
\alias{display.dvh.bands}
\title{Visualizza una banda sul DVH}
\usage{
display.dvh.bands(dvh, plan = NULL, alpha = 0.37, eval.bands = TRUE,
  alpha.color = 0.25, show.plot = TRUE, file.name = NULL, width = 7,
  height = 7, first.is.reference = FALSE, with.mean = TRUE,
  with.median = TRUE)
}
\description{
accetta una lista "dvh.bands" (composta da due dvh: dvh.max
e dvh.min) se first.is.reference=TRUE allora visualizza il
primo dvh della lista come riferimento (non è quindi usato
per calcolare le bande)
}

