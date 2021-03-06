% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reporting-functions.R
\name{generate.report}
\alias{generate.report}
\title{Generate report}
\usage{
generate.report(plan, file.name = "report", N.slice = 6,
  html = FALSE)
}
\arguments{
\item{plan}{the \code{plan} object}

\item{file.name}{the file name}

\item{N.slice}{number of slices to plot for CT and value distributions}

\item{html}{output html file. Note: there is an issue with markdown package that prevents the correct latex equation display in html. Opening the markdown (*.md) file in Rstudio, and then previewing it, overcomes this issue (Rstudio uses an internal markdown to html conversion).}
}
\description{
Generate a report summarizing the plan in a file with the markdown format (it can be converted in html via rstudio or the "markdown" package).
}
