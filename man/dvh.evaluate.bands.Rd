% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dvh-functions.R
\name{dvh.evaluate.bands}
\alias{dvh.evaluate.bands}
\title{Evaluate DVH confidence band}
\usage{
dvh.evaluate.bands(dvh.list, alpha = 0.37)
}
\arguments{
\item{dvh.list}{list of DVHs}

\item{alpha}{the alpha value}
}
\value{
a list of DVHs consisting of:
\item{dvh.mean}{mean}
\item{dvh.median}{median}
\item{dvh.alpha.lo}{lower band}
\item{dvh.alpha.up}{upper band}
\item{dvh.min}{min values}
\item{dvh.max}{max values}
}
\description{
Evaluate the DVH confidence band of a distribution (a list) of DVHs, for a specified alpha (percentile) value.
}
\seealso{
Other DVH: \code{\link{dvh.evaluate.all}},
  \code{\link{dvh.evaluate}}
}
\concept{DVH}
