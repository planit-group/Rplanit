\name{NTCP.LL.LM}
\alias{NTCP.LL.LM}
\title{calcola la log likelihood di complicazioni (NTCP), per un singolo caso, rispetto al dato clinico}
\usage{
NTCP.LL.LM(values, vois, alpha0, beta, m, G = G, s = s, Np = Np, R = R,
  Nf = Nf, voi = voi)
}
\description{
Utilizza il modello lineare (alpha=alpha0 + m*LETd,
beta=cost.). Assume frazioni indipendenti e costanti.
Assume la presenza di matrice di dose fisica e di LETd per
il calcolo LM.
}
\details{
Np -> numero di pazienti R -> numero di pazienti responder
(R=Np*NTCP) Nf -> numero di frazioni
}
\seealso{
Other Likelihood Evaluations: \code{\link{NTCP.LL.cMKM}};
\code{\link{TCP.LL.LM}}; \code{\link{TCP.LL.cMKM}}
}

