\name{tcp.lqp.train}
\alias{tcp.lqp.train}
\title{TCP (Biological Dose, Poissonian Model, Variable fractions)}
\usage{
tcp.lqp.train(dvhs, D50 = 66, R = 10, gam = 2, d = 2, alphaX = NULL,
  betaX = NULL, RBE = 1)
}
\arguments{
  \item{dvhs}{a list of biological dose DVHs}

  \item{D50,R,gam}{tissue specific biological parameters:
  tolerance dose at 50\%, R=alphaX/betaX, clonogenic
  density}

  \item{d}{dose per fraction}

  \item{Nf}{number of fractions}

  \item{RBE}{Relative Biological Effectivenes}
}
\value{
The TCP value
}
\description{
TCP evaluation from the biological dose, using arbitrary
set of fractions.
}
\details{
The total treatment correspond to a "train" of acute doses
distribution (with no time correlation, i.e. full-repair is
assumed). The biological dose distribution, for each voxel,
can be different among fractions. Hence a list of DVH
objects is needed. Each element of the list correspond to a
DHV evaluated for a specific fraction.

By default it is assumed that the \code{dvh.dose} contains
the biological dose dvh [Gy(E)]. Alternatively it could
contain the "physical" dose [Gy] and the biological dose is
evaluated by specifying the RBE per fraction (assumed to be
constant over the irradiated volume).

If alphaX and betaX are not explicity defined, they were
deduced from \code{D50}, \code{gam} and \code{R}.
}
\seealso{
Other General TCP/NTCP Models: \code{\link{ntcp.S}};
\code{\link{ntcp.kallman.train}};
\code{\link{ntcp.kallman}}; \code{\link{tcp.S}};
\code{\link{tcp.lqp}}
}

