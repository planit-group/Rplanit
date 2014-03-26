\name{write.analyze}
\alias{write.analyze}
\title{Write 3D array (analyze)}
\usage{
write.analyze(values = values, variable = variable, file.name = file.name,
  gzipped = FALSE)
}
\arguments{
  \item{values}{the \code{values} object}

  \item{variable}{the variable name}

  \item{file.name}{the file name. Two files will be
  produced: \code{<file.name>.hdr} and
  \code{<file.name>.img}}

  \item{gzipped}{compress the file (boolean, optional)}
}
\description{
Write a \code{values} object in two files using the analyze
format (*.hdr, *.img)
}
\details{
Note: with the analyze format only a single variable array
from \code{values} can be saved in the same files.
}

