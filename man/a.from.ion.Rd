% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/particles.R
\name{a.from.ion}
\alias{a.from.ion}
\title{A from Ion}
\usage{
a.from.ion(ions, with.atomic.mass = FALSE)
}
\arguments{
\item{ions}{A character vector containing the ion names.}

\item{with.atomic.mass}{Uses ion names including the atomic mass number (default FALSE.)}
}
\value{
A numeric vector containing the atomic mass numbers.
}
\description{
Returns the atomic mass number (A) of the atomic nucleus from the name of the atom/ion.
The ion can be specified from its element symbol, optionally including
its atomic mass number (e.g.: "C" or "12C" for carbon ion.)
Note: the available range is from 1H to 20 Ne.
}
\seealso{
Other Particles: \code{\link{Zeff}}, \code{\link{beta.c}},
  \code{\link{ion.from.z}},
  \code{\link{remove.atomic.mass}},
  \code{\link{rest.energy}}, \code{\link{z.from.ion}}
}
\concept{Particles}
