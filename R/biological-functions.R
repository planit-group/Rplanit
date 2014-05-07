# funzioni per calcoli "biologici"...


# VARIABILI DI "SISTEMA" -------------------------------------------------------

#' srim.limits
srim.let.max <- c(H=8.431E-02, He=2.369E-01, Li=4.040E-01, Be=6.052E-01, B=7.825E-01, C=9.757E-01, N=1.159E+00, O=1.376E+00, F=1.609e+00, Ne=1.689e+00)*1000

#' srim.limits
srim.let.min <- c(H=2.060E-04, He=8.243E-04, Li=1.900E-03, Be=3.351E-03, B=5.484E-03, C=8.126E-03, N=1.174E-02, O=1.606E-02, F=1.964e-02, Ne=2.458e-02)*1000

#' srim.limits
srim.es.max <- c(H=10.00E+03, He=2.500e+03, Li=1.428e+03, Be=1.111e+03, B=9.090e+02, C=9.757E-01, N=7.142e+02, O=6.2500e+02, F=6.842680e+02, Ne=6.442118e+02)

#' srim.limits
srim.es.min <- c(H=10.00E-03, He=2.500e-03, Li=1.428e-03, Be=1.111e-03, B=9.090e-04, C=8.333e-04, N=7.142e-04, O=6.2500e-04, F=5.263600e-04, Ne=4.955475e-04)


#' MKM lower limits
mkm.lower <- c(0.1,0.02,1,0.1)

#' MKM upper limits
mkm.upper <- c(1,1,10,1)


# BASIC BIOLOGICAL MODELS ------------------------------------------------------


#' Evaluation of Linear-Quadratic (LQ) Poisson Model parameters
#'
#' It translates R, D50 and gamma parameters to alphaX and betaX.
#' 
#' Note: d.ref è la dose per frazione di riferimento
#' in cui sono fittati i dati per D50, gam.
#' Per il momento d.ref non è usato: si assume una condizione di riferimento
#' standard per D50 e gam (forse... )
#'
#' @param R  the ratio alphaX/betaX (Gy).
#' @param D50  the tolerance dose at 50\%.
#' @param gam  the gamma parameter (associated to the clonogenic density).
#' @param d.ref  the reference dose at which D50 and gam are fitted (usually 2 Gy).
#' @param d  the actual level of dose for the evaluation of alphaX, betaX.
#' 
#' @return A list consisting of:
#' \item{alphaX}{the LQ alpha parameter}
#' \item{betaX}{the LQ beta parameter}
#' 
#' @family Basic Biological Models
#' @export
alpha.beta.lqp <- function(R=R, D50=D50, gam=gam, d.ref=2, d=2)
{
  e <- exp(1)
  D50 <- D50*(1+d.ref/R)/(1+d/R)
  alphaX <- (e*gam-log(log(2)))/(D50*(1+d/R))
  betaX  <- (e*gam-log(log(2)))/(D50*(R+d))
  return(list(alphaX=alphaX, betaX=betaX)) 
}


#' Evaluation of dose-depending RBE
#' 
#' Evaluate the Relative Biological effectiveness (RBE) of the actual radiation.
#' It uses the parameters: alpha, beta, alphaX and betaX
#' 
#' @param alpha,beta LQ parameters for the actual radiation.
#' @param alphaX,betaX LQ parameters for the reference radiation.
#' @param dose dose (Gy)
#' 
#' @return the value of the RBE
#' 
#' @family Basic Biological Models
#' @export
rbe.evaluate <- function(alpha=alpha, beta=beta, dose=dose, alphaX=alphaX, betaX=betaX) {
  R <- alphaX/betaX
  RBE.alpha <- alpha/alphaX
  RBE.beta2 <- beta/betaX
  RBE <- R*( -1 + sqrt( 1 + 4/R * ( RBE.alpha*dose + RBE.beta2*dose^2/R ) ) ) / (2*dose)
  return(RBE)
}


# GENERAL TCP/NTCP MODELS ------------------------------------------------------


#'TCP (Poissonian Model, Biological Dose)
#'
#' Evaluate the Tumor Control Probability (TCP) using the Poissonian Model from the knowledge of the biological dose distribution.
#' 
#' The default parameters are taken from the example in [in Bentzen1997]
#' plus the R=alphaX/betaX parameter for tumoral tissues.
#'  
#' By default it is assumed that the \code{dvh.dose} contains the biological dose dvh [Gy(E)]. Alternatively it could contain the "physical" dose [Gy] and the biological dose is evaluated by specifying the RBE per fraction (assumed to be constant over the irradiated volume).
#'
#' If alphaX and betaX are not explicity defined, they were deduced from \code{D50}, \code{gam} and \code{R}
#' 
#' @param dvh.dose a dvh object in which is stored the dose (biological dose) for the specific volume.
#' @param D50,R,gam tissue specific biological parameters: tolerance dose at 50\%, R=alphaX/betaX, clonogenic density.
#' @param alphaX,betaX alternatively to (D50, R, gam), it is possible to use directly the LQ parameters as biological parameters. These parameters can be vectors.
#' @param d dose per fraction.
#' @param Nf number of fractions.
#' @param G Lea-Catcheside dose-protraction factor (if not define it is assumed complete repair between fractions, G=1/Nf).
#' @param RBE Relative Biological Effectivenes.
#' 
#' @return The TCP value
#' 
#' @family General TCP/NTCP Models
#' @export
#' 
tcp.lqp <- function(dvh.dose=dvh.dose, D50=66, R=10,  gam=2, d=2, alphaX=NULL, betaX=NULL, Nf=30, G=NULL, RBE=1)
{
  # costanti
  e <- exp(1)
 
  # ricava i valori di alphaX e betaX
  if((is.null(alphaX) | is.null(betaX))) {
    alphaX <- alpha.beta.lqp(R, D50, gam, d)$alphaX
    betaX <- alpha.beta.lqp(R, D50, gam, d)$betaX
  }
  #cat('LQ parameters:\n')
  #cat('alphaX =', mean(alphaX), 'Gy^(-1)\n')
  #cat('betaX =', mean(betaX), 'Gy^(-2)\n')

  # Lea-Catcheside dose-protraction factor (assume repair completo tra frazioni)
  if(is.null(G)) {G <- 1/Nf}

  # volume elemetare (frazione):
  dv <- 1/dvh.dose$Nvoxel

  # dose assorbita per frazione
  D.tot <- dvh.dose$value * Nf * RBE  # dose totale assorbita
  #print(length(alphaX))

  S <- exp(-alphaX*D.tot - G*betaX*D.tot*D.tot)
  return( exp(-exp(e*gam)*dv*sum(S)) )
}


#' NTCP (Biological Dose, Relative Seriality Model [Kallman1992a])
#' 
#' It Uses as default parameters those of the "spinal cord" [from Su2010].
#'
#' @param dvh.dose a dvh object in which is stored the dose (biological dose) for the specific volume.
#' @param D50,R,gam tissue specific biological parameters: tolerance dose at 50\%, R=alphaX/betaX, clonogenic density.
#' @param alphaX,betaX alternatively to (D50, R, gam), it is possible to use directly the LQ parameters as biological parameters. These parameters can be vectors.
#' @param s relative seriality parameter of the Kallman's model.
#' @param d dose per fraction.
#' @param Nf number of fractions.
#' @param G Lea-Catcheside dose-protraction factor (if not define it is assumed complete repair between fractions, G=1/Nf).
#' @param RBE Relative Biological Effectivenes.
#' 
#' @return The NTCP value
#' 
#' @family General TCP/NTCP Models
#' @export
#' 
ntcp.kallman <- function(dvh.dose, D50=57, R=3,  gam=6.7, s=1, d=2, alphaX=NULL, betaX=NULL, Nf=30, G=NULL, RBE=1)
{
  # costanti
  e <- exp(1)

  # ricava i valori di alphaX e betaX
  if((is.null(alphaX) | is.null(betaX))) {
    alphaX <- alpha.beta.lqp(R, D50, gam, d)$alphaX
    betaX <- alpha.beta.lqp(R, D50, gam, d)$betaX
  }
  #cat('LQ parameters:\n')
  #cat('alphaX =', mean(alphaX), 'Gy^(-1)\n')
  #cat('betaX =', mean(betaX), 'Gy^(-2)\n')

  # Lea-Catcheside dose-protraction factor (assume repair completo tra frazioni)
  if(is.null(G)) {G <- 1/Nf}

  # volume elemetare (frazione):
  dv <- 1/dvh.dose$Nvoxel

  # dose assorbita per frazione
  D.tot <- dvh.dose$value * Nf * RBE  # dose totale assorbita

  S <- exp(-alphaX*D.tot - G*betaX*D.tot*D.tot)
  return( (1 - prod( (1 - exp(-s*exp(e*gam)*S))^dv ) )^(1/s) )
}


#' TCP (Biological Dose, Poissonian Model, Variable fractions)
#' 
#' TCP evaluation from the biological dose, using arbitrary set of fractions.
#' 
#' The total treatment correspond to a "train" of acute doses distribution (with no time correlation, i.e. full-repair is assumed).
#' The biological dose distribution, for each voxel, can be different among fractions.
#' Hence a list of DVH objects is needed. Each element of the list correspond to a DHV evaluated for a specific fraction.
#' 
#' By default it is assumed that the \code{dvh.dose} contains the biological dose dvh [Gy(E)]. Alternatively it could contain the "physical" dose [Gy] and the biological dose is evaluated by specifying the RBE per fraction (assumed to be constant over the irradiated volume).
#'
#' If alphaX and betaX are not explicity defined, they were deduced from \code{D50}, \code{gam} and \code{R}.
#' 
#' @param dvhs a list of biological dose DVHs
#' @param D50,R,gam tissue specific biological parameters: tolerance dose at 50\%, R=alphaX/betaX, clonogenic density
#' @param d dose per fraction
#' @param Nf number of fractions
#' @param RBE Relative Biological Effectivenes
#' 
#' @return The TCP value
#' 
#' @family General TCP/NTCP Models
#' @export
#' 
tcp.lqp.train <- function(dvhs, D50=66, R=10,  gam=2, d=2, alphaX=NULL, betaX=NULL, RBE=1)
{
  # costanti
  e <- exp(1)
  
  # ricava i valori di alphaX e betaX
  if((is.null(alphaX) | is.null(betaX))) {
    alphaX <- alpha.beta.lqp(R, D50, gam, d)$alphaX
    betaX <- alpha.beta.lqp(R, D50, gam, d)$betaX
  }
  #cat('LQ parameters:\n')
  #cat('alphaX =', mean(alphaX), 'Gy^(-1)\n')
  #cat('betaX =', mean(betaX), 'Gy^(-2)\n')
  
  # numero di frazioni
  Nf <- length(dvhs)
  #message('Number of fractions: ', Nf)
  
  # volume elemetare (frazione). Assume che le divere frazioni abbiano la stessa griglia di voxel
  dv <- 1/dvhs[[1]]$Nvoxel
  
  # sopravvivenza finale per voxel
  S <- dvhs[[1]]$value*0 + 1 # vettore di sopravvivenze finali per ogni voxel
  for(i in 1:Nf) {
    df <- dvhs[[i]]$value * RBE
    s <- exp(-alphaX*df - betaX*df*df)
    S <- S*s # assume scorrelazione temporale tra le frazioni
  }
  
  return( exp(-exp(e*gam)*dv*sum(S)) )
}


#' NTCP (Biological Dose, Relative Seriality Model [Kallman1992a], Variable fractions)
#' 
#' NTCP evaluation from the biological dose, using arbitrary set of fractions.
#' 
#' The total treatment correspond to a "train" of acute doses distribution (with no time correlation, i.e. full-repair is assumed).
#' The biological dose distribution, for each voxel, can be different among fractions.
#' Hence a list of DVH objects is needed. Each element of the list correspond to a DHV evaluated for a specific fraction.
#' 
#' By default it is assumed that the \code{dvh.dose} contains the biological dose dvh [Gy(E)]. Alternatively it could contain the "physical" dose [Gy] and the biological dose is evaluated by specifying the RBE per fraction (assumed to be constant over the irradiated volume).
#'
#' If alphaX and betaX are not explicity defined, they were deduced from \code{D50}, \code{gam} and \code{R}.
#' 
#' @param dvhs a list of biological dose DVHs
#' @param D50,R,gam tissue specific biological parameters: tolerance dose at 50\%, R=alphaX/betaX, clonogenic density
#' @param s relative seriality parameter
#' @param d dose per fraction
#' @param Nf number of fractions
#' @param RBE Relative Biological Effectivenes
#' 
#' @return The NTCP value
#' 
#' @family General TCP/NTCP Models
#' @export
ntcp.kallman.train <- function(dvhs, D50=57, R=3,  gam=6.7, s=1, d=2, alphaX=NULL, betaX=NULL, RBE=1)
{
  # costanti
  e <- exp(1)
  
  # ricava i valori di alphaX e betaX
  if((is.null(alphaX) | is.null(betaX))) {
    alphaX <- alpha.beta.lqp(R, D50, gam, d)$alphaX
    betaX <- alpha.beta.lqp(R, D50, gam, d)$betaX
  }
  #cat('LQ parameters:\n')
  #cat('alphaX =', mean(alphaX), 'Gy^(-1)\n')
  #cat('betaX =', mean(betaX), 'Gy^(-2)\n')

  # numero di frazioni
  Nf <- length(dvhs)
  #message('Number of fractions: ', Nf)
  
  # volume elemetare (frazione). Assume che le divere frazioni abbiano la stessa griglia di voxel
  dv <- 1/dvhs[[1]]$Nvoxel
  
  # sopravvivenza finale per voxel
  S <- dvhs[[1]]$value*0 + 1 # vettore di sopravvivenze finali per ogni voxel
  for(i in 1:Nf) {
    df <- dvhs[[i]]$value * RBE
    sf <- exp(-alphaX*df - betaX*df*df)
    S <- S*sf # assume scorrelazione temporale tra le frazioni
  }
  
  return( (1 - prod( (1 - exp(-s*exp(e*gam)*S))^dv ) )^(1/s) )
}


#' TCP (Poissonian Model)
#'
#' Evaluate TCP directly from the DVH of the "survivals" (for a given endpoint and VOI), without using the biological dose and/or indirect radiobiological parameters for a reference radiation.
#' 
#' It assumes identical and uncorrelated fractions (complete-repair)
#' 
#' @param dvh.survival DVH of the cell "survivals" for a single fraction
#' @param gam radiobiological parameter, related to the clonogenic density
#' @param Nf number of fractions
#' 
#' @return the TCP value.
#' 
#' @family General TCP/NTCP Models
#' @export
tcp.S <- function(dvh.survival=dvh.survival, gam=2, Nf=30) {
  
  # costanti
  e <- exp(1)

  # volume elementare (frazione). Assume che le diverse frazioni abbiano la stessa griglia di voxel
  dv <- 1/dvh.survival$Nvoxel
  
  # Sopravvivenza complessiva (assume completa scorrelazione temporale e che le diverse frazioni siano identiche)
  S <- dvh.survival$value^Nf
  
  return( exp(-exp(e*gam)*dv*sum(S)) )
}


#' NTCP (Relative Seriality Model [Kallman1992a])
#'
#' Evaluate the NTCP directly from the DVH of the "survivals" (for a given endpoint and VOI), without using the biological dose and/or indirect radiobiological parameters for a reference radiation.
#' 
#' It assumes identical and uncorrelated fractions (complete-repair)
#' 
#' @param dvh.survival DVH of the cell "survivals" for a single fraction
#' @param gam radiobiological parameter, related to the clonogenic density
#' @param s relative serialiy parameter
#' @param Nf number of fractions
#' 
#' @return The NTCP value.
#' 
#' @family General TCP/NTCP Models
#' @export
#' 
ntcp.S <- function(dvh.survival, gam=2, s=1, Nf=30) {
  
  # costanti
  e <- exp(1)

  # volume elementare (frazione). Assume che le diverse frazioni abbiano la stessa griglia di voxel
  dv <- 1/dvh.survival$Nvoxel
  
  # Sopravvivenza complessiva (assume completa scorrelazione temporale e che le diverse frazioni siano identiche)
  S <- dvh.survival$value^Nf
  
  return( (1 - prod( (1 - exp(-s*exp(e*gam)*S))^dv ) )^(1/s) )
}


#LEM/MKM MODELS ----------------------------------------------------------------

#' Evaluation of alpha and beta LQ parameter (MKM)
#' 
#' The evaluation is performed for a monoenergetic ion. It uses a C++ MKM rapid implementation for the evaluation ("Survival"). It accepts a single set of MKM parameter (\code{alphaX}, \code{betaX}, \code{rN}, \code{rd}), or alternatively a full range of variability (min,max) for each parameter.
#' 
#' @param alphaX,betaX,rN,rd MKM "biological" parameter associated to a specific biological tissue/cell line:
#' \itemize{
#'   \item{alphaX} LQ alpha parameter of the reference radiation
#'   \item{betaX} LQ beta parameter of the reference radiation
#'   \item{rN} cell nucleus radius [um]
#'   \item{rd} domain radius [um]
#' }
#' 
#' @param alphaX.min,alphaX.max range of variability for parameter \code{alphaX}
#' @param alphaX.N number of step for \code{alphaX}
#' @param betaaX.min,betaX.max range of variability for parameter \code{betaX}
#' @param betaX.N number of step for \code{betaX}
#' @param rN.min,rN.max range of variability for parameter \code{rN}
#' @param rN.N number of step for \code{rN}
#' @param rd.min,rN.max range of variability for parameter \code{rd}
#' @param rd.N number of step for \code{rd}
#' @param cellType name of the tissue/cell line (optional)
#' @param particel type. Available ions: 'H', 'He', 'Li', 'Be', 'B,', 'C', 'N', 'O', 'F', 'Ne'.
#' @param energies vector of energies for the particle [MeV]
#' @param lets vector of LETs for the particle [keV/um]. It is used if \code{energies} is \code{NULL}.
#' 
#' @return a data.frame containing all the information specified including the alpha and beta MKM evaluation (note, in the MKM implementation beta = betaX).
#'
#' @family LEM/MKM Models
#' @export
# utilizza il main main_alpha_beta_parameter._study.cc
# crea una stringa di argomenti da passare a alpha_beta_parameter_study, predisposta
# per il calcolo MKM "rapido".
# In output ritorna un data.frame con i parametri corrispondenti. I valori di
# output vengono passati attraverso il file temporaneo salvato dall'eseguibile, che
# viene cancellato subito dopo. Il file temporaneo è identificato dal nome della
# cellType. Nota: potrebbe essere possibile creare automaticamente N cellType
# e quindi lanciare in parallelo N processi di calcolo...?
alpha.beta.mkm <- function(alphaX=0.1295, betaX=0.03085, rN=4, rd=0.31,
                           alphaX.min=NULL, alphaX.max=NULL, alphaX.N=NULL,
                           betaX.min=NULL, betaX.max=NULL, betaX.N=NULL,
                           rN.min=NULL, rN.max=NULL, rN.N=NULL,
                           rd.min=NULL, rd.max=NULL, rd.N=NULL,
                           cellType=NULL,
                           particleType='H',
                           energies=NULL, lets=NULL,
                           ignore.stdout=TRUE, ignore.stderr=TRUE)
{
  model='MKM'
  calculusType='rapidMKM'
  
  #lem.setenv=get('lem.setenv', envir=dektoolsEnv)
  
  # nome cellType
  if(is.null(cellType)) {cellType <- paste('R', sprintf("%06d", round(runif(1,min=0,max=1e6))), sep='')}

  # valori dei parametri MKM
  if (!is.null(alphaX)) {a <- paste(alphaX, alphaX, 1)} else {a <- paste(alphaX.min, alphaX.max, alphaX.N)}
  if (!is.null(betaX)) {b <- paste(betaX, betaX, 1)} else {b <- paste(betaX.min, betaX.max, betaX.N)}
  if (!is.null(alphaX)) {n <- paste(rN, rN, 1)} else {n <- paste(rN.min, rN.max, rN.N)}
  if (!is.null(alphaX)) {d <- paste(rd, rd, 1)} else {d <- paste(rd.min, rd.max, rd.N)}

  # check per vedere se il calcolo è sulle energie o sul let (e controllo estremi):
  if(!is.null(energies)) {
    energyType <- 'energy'
    e <- paste(energies, collapse=' ')
  } else if(!is.null(lets)) {
    lets <- lets[lets >= srim.let.min[particleType] &  lets <= srim.let.max[particleType]]
    energyType <- 'let'
    e <- paste(lets, collapse=' ')
  } else {
    stop('energies/lets not specified.')
  }

  # output file
  of <- paste(particleType, cellType, model, calculusType, sep='_')
  of <- paste(of, 'csv', sep='.')

  # costruzione linea di comando:
  s.args <- paste(cellType, model, calculusType, a, b, n, d, particleType, energyType, e)
  # cmd <- paste('.', lem.setenv, '; survival_alpha_beta_parameter_study', s.args)
  cmd <- paste('survival_alpha_beta_parameter_study', s.args)
  if(!ignore.stdout) {message(cmd)}

  t <- system.time(system(cmd, ignore.stdout=ignore.stdout, ignore.stderr=ignore.stderr))
  #message('time elapsed: ', t)

  # legge file temporaneo salvato
  out.df <- read.csv(of)

  # cancella file temporaneo
  system(paste('rm', of))

  return(out.df)
}


#' calcolo alpha_beta_parameter_study MKM
#' 
#' utilizza il main main_alpha_beta_parameter._study.cc
#' crea una stringa di argomenti da passare a alpha_beta_parameter_study, predisposta
#' per il calcolo MKM "rapido".
#' In output ritorna la funzione alpha.mkm (funzione interpolazione di R)
#' 
#' @family LEM/MKM Models
#' @export
#' 
alpha.fun.mkm <- function(alphaX=0.1295, betaX=0.03085, rN=4, rd=0.31,
                      cellType=NULL,
                      particleType='H',
                      energies=NULL, lets=NULL,
                      use.limits=FALSE)
{
  model='MKM'
  calculusType='rapidMKM'
  
  #lem.setenv=get('lem.setenv', envir=dektoolsEnv)

  # nome cellType
  if(is.null(cellType)) {
    cellType <- paste('R', sprintf("%06d", round(runif(1,min=0,max=1e6))), sep='')
    #message('cellType = ', cellType)
  }
  
  # valori dei parametri MKM
  a <- paste(alphaX, alphaX, 1)
  b <- paste(betaX, betaX, 1)
  n <- paste(rN, rN, 1)
  d <- paste(rd, rd, 1)

  # check per vedere se il calcolo è sulle energie o sul let:
  if(!is.null(energies)) {
    energyType <- 'energy'
    e <- paste(energies, collapse=' ')
    ee <- energies
  } else {
    lets <- lets[lets >= srim.let.min[particleType] &  lets <= srim.let.max[particleType]]
    energyType <- 'let';
    e <- paste(lets, collapse=' ')
    ee <- lets
  }

  # check per vedere se siamo fuori limiti (restituisce valori alti di alpha, da usare nelle ottimizzazioni)
  if(use.limits) {
    #message('using limits')
    if(!prod(c(alphaX, betaX, rN, rd)>mkm.lower & c(alphaX, betaX, rN, rd)<mkm.upper)) {
      mkm.fun <- approxfun(ee, rep(1e8, length(ee)), method='linear')
      #message(alphaX, betaX, rN, rd, ' out of limits!', sep=' ')
      return(mkm.fun)
    }
  }

  # output file
  of <- paste(particleType, cellType, model, calculusType, sep='_')
  of <- paste(of, 'csv', sep='.')

  # costruzione linea di comando:
  s.args <- paste(cellType, model, calculusType, a, b, n, d, particleType, energyType, e)
  #cmd <- paste('.', lem.setenv, '; survival_alpha_beta_parameter_study', s.args)
  cmd <- paste('survival_alpha_beta_parameter_study', s.args)
  #print(cmd)

  t <- system.time(system(cmd, ignore.stdout=TRUE, ignore.stderr=TRUE))
  #message('time elapsed: ', t)

  # legge file temporaneo salvato
  out.df <- read.csv(of)

  # cancella file temporaneo
  system(paste('rm', of))

  # crea funzione interpolazione
  mkm.fun <- approxfun(ee, out.df$alpha, method='linear')

  #print(c(alphaX, betaX, rN, rd))
  return(mkm.fun)
}


# VALUES BIOLOGICAL EVALUATIONS ------------------------------------------------


#' Calcola la sopravvivenza secondo il modello LM e la aggiunge in values
#' 
#' m, c e beta sono i parametri del modello
#' m -> coefficiente angolare
#' c -> costante
#' beta -> valore di beta (indipendente dal LETd)
#' Nota: l'oggetto values deve avere le due variabili 'Dose[Gy]' e 'DoseAveragedLET[keV/um]'
#' 
#' @family Values Bioogical Evaluations
#' @export
#' 
values.add.survival.LM <- function(values, alpha0, m, beta, add.alpha=TRUE)
{
  
  # identifica gli indici per 'Dose[Gy]' e 'DoseAveragedLET[keV/um]'
  index.dose<- which(values$variables=='Dose[Gy]')
  index.LETd <- which(values$variables=='DoseAveragedLET[keV/um]') 
  # estrae dose e LETd da values
  dose <- values$values[index.dose,,,]
  LETd <- values$values[index.LETd,,,]
    
  # calcola un nuovo array 'Survival.LM' da 'Dose[Gy]' e 'DoseAveragedLET[keV/um]' usando come parametri m, c, e beta
  alpha <- m*LETd + alpha0
  new.array <- exp(-alpha*dose - beta*dose^2)
  
  # aggiungi l'array calcolato all'oggetto values
  values <- add.array.values(values, new.array, variable='Survival.LM')
  
  # aggiungi opzionalmente anche l'array alpha
  if(add.alpha) {values <- add.array.values(values, alpha, variable='Alpha.LM[Gy^(-1)]')}
    
  # ritorna values
  return(values)
}


#' Calcola l'RBE.alpha e l'aggiunge in values
#' 
#' viene passato come argomento il modello:
#' (NULL -> usa "Alpha[Gy^(-1)]")
#' (LM -> usa "Alpha.LM[Gy^(-1)]")
#' (cMKM -> usa "Alpha.cMKM[Gy(-1)]")
#' 
#' assume che la variabile alpha sia già presente in values, atrimenti la calcola
#' esplicitamente.
#' 
#' nota: nell'oggetto values NON deve essere già presente la variabile "Survival"
#' 
#' nota: nei parametri del modello deve essere esplcitata obbligatoriamente anche il beta
#' anche se non viene esplicitamente usatato, viene usato per calcolare le sopravvivenze
#' 
#' @family Values Bioogical Evaluations
#' @export
#'
values.add.RBE.alpha <- function(values, model=NULL, alphaX=0.2, model.LM=NULL, model.cMKM=NULL)
{
  if(is.null(model)) {var.alpha <- 'Alpha[Gy^(-1)]'; var.rbe <- 'RBE.alpha'}
  else if(model=='LM') {var.alpha <- 'Alpha.LM[Gy^(-1)]'; var.rbe <- 'RBE.alpha.LM'}
  else if(model=='cMKM') {var.alpha <- 'Alpha.cMKM[Gy^(-1)]'; var.rbe <- 'RBE.alpha.cMKM'}
  v.alpha <- which(values$variables==var.alpha)

  # calcola alpha ex novo
  if(length(v.alpha)==0) {
    if(model=='LM') {
      message('Evaluating RBE.alpha using Model: LM')
      values <- values.add.survival.LM(values, alpha0=model.LM$alpha0, m=model.LM$m, beta=model.LM$beta)
    } else if(model=='cMKM') {
      message('Evaluating RBE.alpha using Model: LM')
      values <- values.add.survivals.LM(values,
                                        alphaX=model.cMKM$alphaX, betaX=model.cMKM$betaX,
                                        rN=model.cMKM$rN, rd=model.cMKM$rd,
                                        beta=model.cMKM$beta,
                                        particleType=as.charcater(model.cMKM$particleType),
                                        beta=beta, add.alpha=TRUE)
    }
  }
  v.alpha <- which(values$variables==var.alpha)

  rbe.array <- values$values[v.alpha,,,]/alphaX
  values <- add.array.values(values=values, new.array=rbe.array, variable=var.rbe)

  return(values)
}


#' Calcola la sopravvivenza secondo il modello cMKM e la aggiunge in values
#' 
#' alphaX, betaX, rN, rd sono i parametri del modello MKM
#' particleType -> ione di riferimento
#' beta -> valore di beta (indipendente dal LETd)
#' Nota: l'oggetto values deve avere le due variabili 'Dose[Gy]' e 'DoseAveragedLET[keV/um]'
#' 
#' @family Values Bioogical Evaluations
#' @export
#'
values.add.survival.cMKM <- function(values=values, alphaX=alphaX, betaX=betaX, rN=rN, rd=rd, particleType=particleType, beta=beta, add.alpha=TRUE)
{
  
  # identifica gli indici per 'Dose[Gy]' e 'DoseAveragedLET[keV/um]'
  index.dose<- which(values$variables=='Dose[Gy]')
  index.LETd <- which(values$variables=='DoseAveragedLET[keV/um]') 
  # estrae dose e LETd da values
  dose <- values$values[index.dose,,,]
  LETd <- values$values[index.LETd,,,]

  # vettore di lets da interpolare
  lets <- seq(min(LETd, na.rm=TRUE), max(LETd, na.rm=TRUE), length.out=30)
   
  # calcola un nuovo array 'Survival.cMKM' da 'Dose[Gy]' e 'DoseAveragedLET[keV/um]'
  alpha.mkm <- alpha.fun.mkm(alphaX=alphaX, betaX=betaX, rN=rN, rd=rd,
                             particleType=particleType,
                             lets=lets,
                             use.limits=FALSE)

  alpha <- alpha.mkm(LETd)
  new.array <- exp(-alpha*dose - beta*dose^2)
  
  # aggiungi l'array calcolato all'oggetto values
  values <- add.array.values(values=values, new.array=new.array, variable='Survival.cMKM')
  
  # aggiungi opzionalmente anche l'array alpha
  if(add.alpha) {values <- add.array.values(values, alpha, variable='Alpha.cMKM[Gy^(-1)]')}
    
  # ritorna values
  return(values)
}


#' calcola l'RBE e l'agiunge in values a partire da alpha e beta e dose
#' 
#' @family Values Bioogical Evaluations
#' @export
#'
values.add.RBE <- function(values, model=NULL, alphaX=0.2, betaX=0.02, model.LM=NULL, model.cMKM=NULL, biological.dose=FALSE)
{
  if(is.null(model)) {var.alpha <- 'Alpha[Gy^(-1)]'; var.beta <- 'Beta[Gy^(-2)]'; var.rbe <- 'RBE'}
  else if(model=='LM') {var.alpha <- 'Alpha.LM[Gy^(-1)]'; var.beta <- 'Beta.LM[Gy^(-2)]'; var.rbe <- 'RBE.LM'}
    else if(model=='cMKM') {var.alpha <- 'Alpha.cMKM[Gy^(-1)]'; var.beta <- 'Beta.cMKM[Gy^(-2)]'; var.rbe <- 'RBE.cMKM'}

    v.alpha <- which(values$variables==var.alpha)
    v.beta <- which(values$variables==var.beta)

    # calcola alpha e beta ex novo, e quindi l'rbe
    if(length(v.alpha)+length(v.beta)==0) {
      # modello LM
      if(model=='LM') {
        message('Evaluating RBE using Model: LM')
        index.dose <- which(values$variables=='Dose[Gy]')
        index.LETd <- which(values$variables=='DoseAveragedLET[keV/um]')
	      LETd <- values$values[index.LETd,,,]
	      dose <- values$values[index.dose,,,]
	      alpha <- model.LM$m*LETd + model.LM$alpha0
        beta <- alpha*0 + model.LM$beta
	      rbe <- rbe.evaluate(alpha=alpha, beta=beta, dose=dose, alphaX=alphaX, betaX=betaX)
	      values <- add.array.values(values=values, new.array=rbe, variable='RBE.LM')
        if(biological.dose) {
          dbio <- dose*rbe
          values <- add.array.values(values=values, new.array=dbio, variable='BiologicalDose.LM[Gy(RBE)]')
        }
      }
      # modello cMKM
      else if(model=='cMKM') {
        message('Evaluating RBE using Model: cMKM')
        index.dose <- which(values$variables=='Dose[Gy]')
        index.LETd <- which(values$variables=='DoseAveragedLET[keV/um]')
	      LETd <- values$values[index.LETd,,,]
	      dose <- values$values[index.dose,,,]
	      lets <- seq(min(LETd, na.rm=TRUE), max(LETd, na.rm=TRUE), length.out=30)
	      alpha.mkm <- alpha.fun.mkm(alphaX=model.cMKM$alphaX, betaX=model.cMKM$betaX, rN=model.cMKM$rN, rd=model.cMKM$rd,
				                           cellType='RCell', particleType=as.character(model.cMKM$particleType), lets=lets, use.limits=FALSE)
	      alpha <- alpha.mkm(LETd)
	      beta <- alpha*0 + model.cMKM$beta
	      rbe <- rbe.evaluate(alpha=alpha, beta=beta, dose=dose, alphaX=alphaX, betaX=betaX)
	      values <- add.array.values(values=values, new.array=rbe, variable='RBE.cMKM')
        if(biological.dose) {
          dbio <- dose*rbe
          values <- add.array.values(values=values, new.array=dbio, variable='BiologicalDose.cMKM[Gy(RBE)]')
        }
      }
    }

    return(values)
}


# LM/cMKM TCP/NTCP MODELS ------------------------------------------------------


#' calcola TCP con modello LM
#' 
#' @family LM/cMKM TCP/NTCP Models
#' @export
#'
TCP.LM <- function(values, vois, voi='PTV', alpha0, beta, m, G=G, Nf=Nf){
  values.S <- values.add.survival.LM(values=values, alpha0=alpha0, m=m, beta=beta)
  dvh.surv <- dvh.evaluate(values=values.S, variable='Survival.LM', voi=voi, vois=vois)
  TCP <- tcp.S(dvh.survival=dvh.surv, gam=G, Nf=Nf)
  return(TCP)
}


#' calcola TCP con modello cMKM
#' @family LM/cMKM TCP/NTCP Models
#' @export
#'
TCP.cMKM <- function(values=values, vois=vois, voi='PTV', alphaX=alphaX, betaX=betaX, rN=rN, rd=rd, beta=beta, particleType=particleType, G=G, Nf=Nf){
  values.S <- values.add.survival.cMKM(values=values, alphaX=alphaX, betaX=betaX, rN=rN, rd=rd, beta=beta, particleType=particleType)
  dvh.surv <- dvh.evaluate(values=values.S, variable='Survival.cMKM', voi=voi, vois=vois)
  TCP <- tcp.S(dvh.survival=dvh.surv, gam=G, Nf=Nf)
  return(TCP)
}


#' calcola NTCP con modello LM
#' @family LM/cMKM TCP/NTCP Models
#' @export
#'
NTCP.LM <- function(values, vois, voi='PTV', alpha0, beta, m, G=G, s=s, Nf=Nf){
  values.S <- values.add.survival.LM(values=values, alpha0=alpha0, m=m, beta=beta)
  dvh.surv <- dvh.evaluate(values=values.S, variable='Survival.LM', voi=voi, vois=vois)
  NTCP <- ntcp.S(dvh.survival=dvh.surv, gam=G, s=s, Nf=Nf)
  return(NTCP)
}


#' calcola NTCP con modello cMKM
#' 
#' @family LM/cMKM TCP/NTCP Models
#' @export
#'
NTCP.cMKM <- function(values=values, vois=vois, voi='PTV', alphaX=alphaX, betaX=betaX, rN=rN, rd=rd, beta=beta, particleType=particleType, G=G, s=s, Nf=Nf){
  values.S <- values.add.survival.cMKM(values=values, alphaX=alphaX, betaX=betaX, rN=rN, rd=rd, beta=beta, particleType=particleType)
  dvh.surv <- dvh.evaluate(values=values.S, variable='Survival.cMKM', voi=voi, vois=vois)
  NTCP <- ntcp.S(dvh.survival=dvh.surv, gam=G, s=s, Nf=Nf)
  return(NTCP)
}


# LIKELIHOOD EVALUATIONS -------------------------------------------------------


#' calcola la log likelihood di cura (TCP), per un singolo caso, rispetto al dato clinico
#' 
#' Utilizza il modello lineare (alpha=alpha0 + m*LETd, beta=cost.).
#' Assume frazioni indipendenti e costanti.
#' Assume la presenza di matrice di dose fisica e di LETd per il calcolo LM.
#'
#' Np -> numero di pazienti
#' R -> numero di pazienti responder (R=Np*TCP)
#' Nf -> numero di frazioni
#' 
#' @family Likelihood Evaluations
#' @export
#'
TCP.LL.LM <- function(values, vois, voi='PTV', alpha0, beta, m, G=G, Np=Np, R=R, Nf=Nf){
  values.S <- values.add.survival.LM(values=values, alpha0=alpha0, m=m, beta=beta)
  dvh.surv <- dvh.evaluate(values=values.S, variable='Survival.LM', voi=voi, vois=vois)
  TCP <- tcp.S(dvh.survival=dvh.surv, gam=G, Nf=Nf)
  LL <- log(TCP)*R+log(1-TCP)*(Np-R)
  return(LL)
}


#' calcola la log likelihood di cura (TCP), per un singolo caso, rispetto al dato clinico
#' 
#' Utilizza il modello cMKM.
#' Assume frazioni indipendenti e costanti.
#' Assume la presenza di matrice di dose fisica e di LETd per il calcolo LM.
#'
#' Np -> numero di pazienti
#' R -> numero di pazienti responder (R=Np*TCP)
#' Nf -> numero di frazioni
#' @family Likelihood Evaluations
#' @export
#'
TCP.LL.cMKM <- function(values, vois, voi='PTV',
                      alphaX, betaX, rN, rd, beta, particleType,
                      G, Np, R, Nf){
  values.S <- values.add.survival.cMKM(values=values, alphaX=alphaX, betaX=betaX, rN=rN, rd=rd, particleType=particleType, beta=beta)
  dvh.surv <- dvh.evaluate(values=values.S, variable='Survival.LM', voi=voi, vois=vois)
  TCP <- tcp.S(dvh.survival=dvh.surv, gam=G, Nf=Nf)
  LL <- log(TCP)*R+log(1-TCP)*(Np-R)
  return(LL)
}


#' calcola la log likelihood di complicazioni (NTCP), per un singolo caso, rispetto al dato clinico
#' 
#' Utilizza il modello lineare (alpha=alpha0 + m*LETd, beta=cost.).
#' Assume frazioni indipendenti e costanti.
#' Assume la presenza di matrice di dose fisica e di LETd per il calcolo LM.
#'
#' Np -> numero di pazienti
#' R -> numero di pazienti responder (R=Np*NTCP)
#' Nf -> numero di frazioni
#' @family Likelihood Evaluations
#' @export
#'
NTCP.LL.LM <- function(values, vois, alpha0, beta, m, G=G, s=s, Np=Np, R=R, Nf=Nf, voi=voi){
  values.S <- values.add.survival.LM(values=values, alpha0=alpha0, m=m, beta=beta)
  dvh.surv <- dvh.evaluate(values=values.S, variable='Survival.LM', voi=voi, vois=vois)
  NTCP <- ntcp.S(dvh.survival=dvh.surv, gam=G, s=s, Nf=Nf)
  LL <- log(NTCP)*R+log(1-NTCP)*(Np-R)
  if(is.infinite(LL)) {message('NTCP: ', NTCP)}
  return(LL)
}


#' calcola la log likelihood di complicazioni (NTCP), per un singolo caso, rispetto al dato clinico
#' 
#' Utilizza il modello cMKM.
#' Assume frazioni indipendenti e costanti.
#' Assume la presenza di matrice di dose fisica e di LETd per il calcolo LM.
#'
#' Np -> numero di pazienti
#' R -> numero di pazienti responder (R=Np*NTCP)
#' Nf -> numero di frazioni
#' @family Likelihood Evaluations
#' @export
#'
NTCP.LL.cMKM <- function(values, vois,
                         alphaX, betaX, rN, rd, beta, particleType,
                         G=G, s=s, Np=Np, R=R, Nf=Nf, voi=voi){
  values.S <- values.add.survival.cMKM(values=values, alphaX=alphaX, betaX=betaX, rN=rN, rd=rd, particleType=particleType, beta=beta)
  dvh.surv <- dvh.evaluate(values=values.S, variable='Survival.LM', voi=voi, vois=vois)
  NTCP <- ntcp.S(dvh.survival=dvh.surv, gam=G, s=s, Nf=Nf)
  LL <- log(NTCP)*R+log(1-NTCP)*(Np-R)
  if(is.infinite(LL)) {message('NTCP: ', NTCP)}
  return(LL)
} 
                                       
