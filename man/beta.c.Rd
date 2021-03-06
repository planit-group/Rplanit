% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/particles.R
\name{beta.c}
\alias{beta.c}
\title{Evaluate beta (velocity/c)}
\usage{
beta.c(Ec, A)
}
\arguments{
\item{Ec}{A numeric vector containing the kinetic energies (MeV)}

\item{A}{A numeric vector containing the atomic mass number (or mass in atomic mass units).}
}
\value{
the velocity in c units.
}
\description{
Evaluate the velocity in c units from the kinetic energy of the particle and its atomic mass number.
}
\seealso{
Other Particles: \code{\link{Zeff}},
  \code{\link{a.from.ion}}, \code{\link{ion.from.z}},
  \code{\link{remove.atomic.mass}},
  \code{\link{rest.energy}}, \code{\link{z.from.ion}}
}
\concept{Particles}
