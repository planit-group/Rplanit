% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biological-functions.R
\name{alpha.beta.mkm}
\alias{alpha.beta.mkm}
\title{Evaluation of alpha and beta LQ parameter (fast MKM)}
\usage{
alpha.beta.mkm(alphaX = 0.1295, betaX = 0.03085, rN = 4, rd = 0.31,
  alphaX.min = NULL, alphaX.max = NULL, alphaX.N = NULL,
  betaX.min = NULL, betaX.max = NULL, betaX.N = NULL,
  rN.min = NULL, rN.max = NULL, rN.N = NULL, rd.min = NULL,
  rd.max = NULL, rd.N = NULL, cellType = NULL, particleType = "H",
  energies = NULL, lets = NULL, calculusType = "rapidMKM",
  ignore.stdout = TRUE, ignore.stderr = TRUE)
}
\arguments{
\item{alphaX, betaX, rN, rd}{MKM "biological" parameter associated to a specific biological tissue/cell line:
\itemize{
  \item{alphaX} LQ alpha parameter of the reference radiation [Gy^-1]
  \item{betaX} LQ beta parameter of the reference radiation [Gy^-2]
  \item{rN} cell nucleus radius [um]
  \item{rd} domain radius [um]
}}

\item{alphaX.min, alphaX.max}{range of variability for parameter \code{alphaX}}

\item{alphaX.N}{number of step for \code{alphaX}}

\item{betaX.min, betaX.max}{range of variability for parameter \code{betaX}}

\item{betaX.N}{number of step for \code{betaX}}

\item{rN.min, rN.max}{range of variability for parameter \code{rN}}

\item{rN.N}{number of step for \code{rN}}

\item{rd.min, rd.max}{range of variability for parameter \code{rd}}

\item{rd.N}{number of step for \code{rd}}

\item{cellType}{name of the tissue/cell line (optional)}

\item{energies}{vector of energies for the particle [MeV]}

\item{lets}{vector of LETs for the particle [keV/um]. It is used if \code{energies} is \code{NULL}.}

\item{calculusType}{the type of the evaluations options are:
\itemize{
  \item{rapidKase} implementation of MKM as described in (Kase2008)
  \item{rapidMKM} more correct implementation for the MKM from the one described in (Kase2008)
  \item{newMKMKase} "rapidKase" implementation +  non poissionian correction for beta
  \item{newMKM} "rapidMKM" implementation +  non poissionian correction for beta
  \item{slow_alphaIon_betaIon} slow Monte Carlo evaluation (compatible with MKM, LEMI, LEMII and LEMIII).
}}

\item{particel}{type. Available ions: 'H', 'He', 'Li', 'Be', 'B,', 'C', 'N', 'O', 'F', 'Ne'.}
}
\value{
a data.frame containing all the information specified including the alpha and beta MKM evaluation (note, in the MKM implementation beta = betaX).
}
\description{
The evaluation is performed for a monoenergetic ion. It uses a C++ MKM "rapid" implementation for the evaluation. It accepts a single set of MKM parameter (\code{alphaX}, \code{betaX}, \code{rN}, \code{rd}), or alternatively a full range of variability (min,max) for each parameter.
}
\seealso{
Other LEM/MKM Models: \code{\link{alpha.beta.ion.range}},
  \code{\link{alpha.fun.mkm}}
}
\concept{LEM/MKM Models}
