\name{alpha.fun.mkm}
\alias{alpha.fun.mkm}
\title{calcolo alpha_beta_parameter_study MKM}
\usage{
alpha.fun.mkm(alphaX = 0.1295, betaX = 0.03085, rN = 4, rd = 0.31,
  cellType = NULL, particleType = "H", energies = NULL, lets = NULL,
  use.limits = FALSE)
}
\description{
utilizza il main main_alpha_beta_parameter._study.cc crea
una stringa di argomenti da passare a
alpha_beta_parameter_study, predisposta per il calcolo MKM
"rapido". In output ritorna la funzione alpha.mkm (funzione
interpolazione di R)
}
\seealso{
Other LEM/MKM Models: \code{\link{alpha.beta.mkm}}
}

