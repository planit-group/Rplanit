\name{alpha.beta.mkm}
\alias{alpha.beta.mkm}
\title{Evaluation of alpha and beta LQ parameter (MKM)}
\usage{
alpha.beta.mkm(alphaX = 0.1295, betaX = 0.03085, rN = 4, rd = 0.31,
  alphaX.min = NULL, alphaX.max = NULL, alphaX.N = NULL,
  betaX.min = NULL, betaX.max = NULL, betaX.N = NULL, rN.min = NULL,
  rN.max = NULL, rN.N = NULL, rd.min = NULL, rd.max = NULL,
  rd.N = NULL, cellType = NULL, particleType = "H", energies = NULL,
  lets = NULL)
}
\arguments{
  \item{alphaX,betaX,rN,rd}{MKM "biological" parameter
  associated to a specific biological tissue/cell line:
  \itemize{ \item{alphaX} LQ alpha parameter of the
  reference radiation \item{betaX} LQ beta parameter of the
  reference radiation \item{rN} cell nucleus radius [um]
  \item{rd} domain radius [um] }}

  \item{alphaX.min,alphaX.max}{range of variability for
  parameter \code{alphaX}}

  \item{alphaX.N}{number of step for \code{alphaX}}

  \item{betaaX.min,betaX.max}{range of variability for
  parameter \code{betaX}}

  \item{betaX.N}{number of step for \code{betaX}}

  \item{rN.min,rN.max}{range of variability for parameter
  \code{rN}}

  \item{rN.N}{number of step for \code{rN}}

  \item{rd.min,rN.max}{range of variability for parameter
  \code{rd}}

  \item{rd.N}{number of step for \code{rd}}

  \item{cellType}{name of the tissue/cell line (optional)}

  \item{particel}{type. Available ions: 'H', 'He', 'Li',
  'Be', 'B,', 'C', 'N', 'O', 'F', 'Ne'.}

  \item{energies}{vector of energies for the particle
  [MeV]}

  \item{lets}{vector of LETs for the particle [keV/um]. It
  is used if \code{energies} is \code{NULL}.}
}
\value{
a data.frame containing all the information specified
including the alpha and beta MKM evaluation (note, in the
MKM implementation beta = betaX).
}
\description{
The evaluation is performed for a monoenergetic ion. It
uses a C++ MKM rapid implementation for the evaluation
("Survival"). It accepts a single set of MKM parameter
(\code{alphaX}, \code{betaX}, \code{rN}, \code{rd}), or
alternatively a full range of variability (min,max) for
each parameter.
}
\seealso{
Other LEM/MKM Models: \code{\link{alpha.fun.mkm}}
}

