% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biological-functions.R
\name{rbe.evaluate}
\alias{rbe.evaluate}
\title{Evaluation of dose-depending RBE}
\usage{
rbe.evaluate(alpha = alpha, beta = beta, dose = dose,
  alphaX = alphaX, betaX = betaX)
}
\arguments{
\item{alpha, beta}{LQ parameters for the actual radiation.}

\item{dose}{dose (Gy)}

\item{alphaX, betaX}{LQ parameters for the reference radiation.}
}
\value{
the value of the RBE
}
\description{
Evaluate the Relative Biological effectiveness (RBE) of the actual radiation.
It uses the parameters: alpha, beta, alphaX and betaX
}
\seealso{
Other Basic Biological Models: \code{\link{alpha.beta.lqp}}
}
\concept{Basic Biological Models}
