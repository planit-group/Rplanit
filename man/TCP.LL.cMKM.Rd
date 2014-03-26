\name{TCP.LL.cMKM}
\alias{TCP.LL.cMKM}
\title{calcola la log likelihood di cura (TCP), per un singolo caso, rispetto al dato clinico}
\usage{
TCP.LL.cMKM(values, vois, voi = "PTV", alphaX, betaX, rN, rd, beta,
  particleType, G, Np, R, Nf)
}
\description{
Utilizza il modello cMKM. Assume frazioni indipendenti e
costanti. Assume la presenza di matrice di dose fisica e di
LETd per il calcolo LM.
}
\details{
Np -> numero di pazienti R -> numero di pazienti responder
(R=Np*TCP) Nf -> numero di frazioni
}
\seealso{
Other Likelihood Evaluations: \code{\link{NTCP.LL.LM}};
\code{\link{NTCP.LL.cMKM}}; \code{\link{TCP.LL.LM}}
}

