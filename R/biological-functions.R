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

#' Evaluation of alpha and beta LQ parameter (fast MKM)
#'
#' The evaluation is performed for a monoenergetic ion. It uses a C++ MKM "rapid" implementation for the evaluation. It accepts a single set of MKM parameter (\code{alphaX}, \code{betaX}, \code{rN}, \code{rd}), or alternatively a full range of variability (min,max) for each parameter.
#'
#' @param alphaX,betaX,rN,rd MKM "biological" parameter associated to a specific biological tissue/cell line:
#' \itemize{
#'   \item{alphaX} LQ alpha parameter of the reference radiation [Gy^-1]
#'   \item{betaX} LQ beta parameter of the reference radiation [Gy^-2]
#'   \item{rN} cell nucleus radius [um]
#'   \item{rd} domain radius [um]
#' }
#'
#' @param alphaX.min,alphaX.max range of variability for parameter \code{alphaX}
#' @param alphaX.N number of step for \code{alphaX}
#' @param betaX.min,betaX.max range of variability for parameter \code{betaX}
#' @param betaX.N number of step for \code{betaX}
#' @param rN.min,rN.max range of variability for parameter \code{rN}
#' @param rN.N number of step for \code{rN}
#' @param rd.min,rd.max range of variability for parameter \code{rd}
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
# utilizza il main main_alpha_beta_parameter_study.cc
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
  model <- 'MKM'
  calculusType <- 'rapidMKM'
  precision <- 1
  
  #lem.setenv=get('lem.setenv', envir=dektoolsEnv)
  
  # nome cellType
  if(is.null(cellType)) {cellType <- paste('R', sprintf("%06d", round(runif(1,min=0,max=1e6))), sep='')}
  
  # valori dei parametri MKM
  if (!is.null(alphaX)) {a <- paste(alphaX, alphaX, 1)} else {a <- paste(alphaX.min, alphaX.max, alphaX.N)}
  if (!is.null(betaX)) {b <- paste(betaX, betaX, 1)} else {b <- paste(betaX.min, betaX.max, betaX.N)}
  if (!is.null(rN)) {n <- paste(rN, rN, 1)} else {n <- paste(rN.min, rN.max, rN.N)}
  if (!is.null(rd)) {d <- paste(rd, rd, 1)} else {d <- paste(rd.min, rd.max, rd.N)}
  
  # check per vedere se il calcolo è sulle energie o sul let (e controllo estremi):
  if(length(energies)>0) {
    energyType <- 'energy'
    e <- paste(energies, collapse=' ')
  } else if(length(lets)>0) {
    lets <- lets[lets >= srim.let.min[particleType] &  lets <= srim.let.max[particleType]]
    energyType <- 'let'
    e <- paste(lets, collapse=' ')
  } else if(length(energies)==0 & length(lets)==0){
    stop('energies/lets not specified.')
  }
  
  # output file
  of <- paste(particleType, cellType, model, calculusType, sep='_')
  of <- paste(of, 'csv', sep='.')
  
  # costruzione linea di comando:
  s.args <- paste(cellType, model, calculusType, a, b, n, d, particleType, precision, energyType, e)
  # cmd <- paste('.', lem.setenv, '; survival_alpha_beta_parameter_study', s.args)
  cmd <- paste('survival', s.args)
  if(!ignore.stdout) {message(cmd)}
  
  t <- system.time(system(cmd, ignore.stdout=ignore.stdout, ignore.stderr=ignore.stderr))
  #message('time elapsed: ', t)
  
  # legge file temporaneo salvato
  out.df <- read.csv(of)
  
  # cancella file temporaneo
  system(paste('rm', of))
  
  return(out.df)
}

#' Evaluation of ion alpha and beta LQ parameters
#'
#' The evaluation is performed for a monoenergetic ion. It uses a C++ implementation for the evaluation. It accepts a single set of parameter depending on the choosen radiobiological model, or alternatively a full range of variability (min,max) for each parameter.
#'
#' @param model the name of the model (options: 'MKM', 'LEMI', 'LEMII', 'LEMIII')
#' @param model.parameters a dataframe containing the "biological" parameter associated to a specific biological tissue/cell line. In the case of the MKM:
#' \itemize{
#'   \item{alpha0} LQ alpha parameter of the reference radiation [Gy^-1]
#'   \item{beta0} LQ beta parameter of the reference radiation [Gy^-2]
#'   \item{rN} cell nucleus radius [um]
#'   \item{rd} domain radius [um]
#' }
#' In the case of LEMI, LEMII, LEMIII:
#' \itemize{
#'   \item{alpha0} LQ alpha parameter of the reference radiation [Gy^-1]
#'   \item{beta0} LQ beta parameter of the reference radiation [Gy^-2]
#'   \item{rN} cell nucleus radius [um]
#'   \item{Dt} threshold dose [Gy]
#' }
#' Optionally, instead of the specific values, the range of the parameters can be specified in the same data.frame. In the case of MKM:
#' \itemize{
#'   \item{alpha0.min, alpha0.max, alpha0.N} range of variability for parameter \code{alpha0} and number of steps
#'   \item{beta0.min, beta0.max, beta0.N} range of variability for parameter \code{beta0} and number of steps
#'   \item{rN.min, rN.max, rN.N} range of variability for parameter \code{rN} and number of steps
#'   \item{rd.min, rN.max} domain radius [um]
#' }
#'
#' @param cellType name of the tissue/cell line (optional)
#' @param particle available particles: 'H', 'He', 'Li', 'Be', 'B,', 'C', 'N', 'O', 'F', 'Ne'.
#' @param energies vector of energies for the particle [MeV]
#' @param lets vector of LETs for the particle [keV/um]. It is used if \code{energies} is \code{NULL}.
#' @param calculusType the type of the evaluations options are:
#' \itemize{
#'   \item{rapidMKM} fast implementation for the MKM (Kase2008)
#'   \item{rapidScholz} fast implementation of the LEMI, LEMII, LEMIII (Scholz2000)
#'   \item{rapidRusso} fast implementation of the LEMI, LEMII, LEMIII, more accurate than the \code{rapidScholz} (Russo2010)
#'   \item{slow_alphaIon_betaIon} slow Monte Carlo evaluation (compatible with MKM, LEMI, LEMII and LEMIII).
#' }
#' @param precision used only for the Monte Carlo evaluation. If precision < 1.0 it represents the target relative standard deviation to be reached in the evaluations of the cell survivals. If precision >= 1 it represents the number of cells to be simulated.
#' @param get.raw.data If calculusType=slow_alphaIon_betaIon, returns the raw data (survival vs dose for each simulated cell) from which alpha and beta parameters can be extracted, otherwhise it is ignored.
#' 
#' @return a data.frame containing all the information specified including the alpha and beta MKM evaluation (note, in the MKM implementation beta = betaX).
#'
#' @family LEM/MKM Models
#' @export
alpha.beta.ion.range <- function(model='MKM',
                                 model.parameters=data.frame(alpha0=0.1295, beta0=0.03085, rN=4, rd=0.31),
                                 cell.name=NULL,
                                 particle='H',
                                 energies=NULL, lets=NULL,
                                 calculusType='rapidMKM', precision=0.5,
                                 ignore.stdout=TRUE, ignore.stderr=TRUE,
                                 get.raw.data=FALSE, remove.temp.files=TRUE)
{
  # nome cellType
  if(is.null(cell.name)) {cell.name <- paste('R', sprintf("%06d", round(runif(1,min=0,max=1e6))), sep='')}
  
  # nomi dei parametri
  if(model=='MKM') {
    parameter.names <- c('alpha0', 'beta0', 'rN', 'rd')
    parameter.names.min <- c('alpha0.min', 'beta0.min', 'rN.min', 'rd.min')
    parameter.names.max <- c('alpha0.max', 'beta0.max', 'rN.max', 'rd.max')
    parameter.names.N <- c('alpha0.N', 'beta0.N', 'rN.N', 'rd.N')
  }
  else if(model=='LEMI' | model=='LEMII' | model=='LEMIII') {
    parameter.names <- c('alpha0', 'beta0', 'rN', 'Dt')
    parameter.names.min <- c('alpha0.min', 'beta0.min', 'rN.min', 'Dt.min')
    parameter.names.max <- c('alpha0.max', 'beta0.max', 'rN.max', 'Dt.max')
    parameter.names.N <- c('alpha0.N', 'beta0.N', 'rN.N', 'Dt.N')
  } else {
    stop('Model', model, 'not defined.')
  }
  
  # recupera parametri
  if(all(parameter.names %in% names(model.parameters))) {
    p <- model.parameters[parameter.names]
    p.min <- p.max <- p.N <- NULL
  } else if(all(c(parameter.names.min, parameter.names.max, parameter.names.N) %in% names(model.parameters))) {
    p.min <- model.parameters[parameter.names.min];
    p.max <- model.parameters[parameter.names.max];
    p.N <- model.parameters[parameter.names.N];
    p <- NULL
  } else {
    stop('Parameter not well defined. (', names(model.parameters), ')')
  }
  
  # valori dei parametri da passare a pure-survival
  if (!is.null(p[1])) {a <- paste(p[1], p[1], 1)} else {a <- paste(p.min[1], p.max[1], p.N[1])}
  if (!is.null(p[2])) {b <- paste(p[2], p[2], 1)} else {b <- paste(p.min[2], p.max[2], p.N[2])}
  if (!is.null(p[3])) {n <- paste(p[3], p[3], 1)} else {n <- paste(p.min[3], p.max[3], p.N[3])}
  if (!is.null(p[4])) {d <- paste(p[4], p[4], 1)} else {d <- paste(p.min[4], p.max[4], p.N[4])}
  
  # check per vedere se il calcolo è sulle energie o sul let (e controllo estremi):
  if(length(energies)>0) {
    energyType <- 'energy'
    e <- paste(energies, collapse=' ')
  } else if(length(lets)>0) {
    lets <- lets[lets >= srim.let.min[particle] &  lets <= srim.let.max[particle]]
    energyType <- 'let'
    e <- paste(lets, collapse=' ')
  } else if(length(energies)==0 & length(lets)==0){
    stop('energies/lets not specified.')
  }
  
  # output file
  of <- paste(particle, cell.name, model, calculusType, sep='_')
  of <- paste(of, 'csv', sep='.')
  
  # creazione cartella per calcolo MC
  if(calculusType=='slow_alphaIon_betaIon') {dir.create(paste0('survival-data-', particle), showWarnings = FALSE)}
  
  # costruzione linea di comando:
  s.args <- paste(cell.name, model, calculusType, a, b, n, d, particle, precision, energyType, e)
  # cmd <- paste('.', lem.setenv, '; survival_alpha_beta_parameter_study', s.args)
  cmd <- paste('survival', s.args)
  if(!ignore.stdout) {message(cmd)}
  
  t <- system.time(system(cmd, ignore.stdout=ignore.stdout, ignore.stderr=ignore.stderr))
  #message('time elapsed: ', t)
  
  # legge file temporaneo salvato
  out.df <- read.csv(of)
  out.df$cell <- cell.name 
  out.df$particle <- particle
  
  # cancella file temporaneo
  if(remove.temp.files) {
    #system(paste('rm', of))
    unlink(of)
    if(calculusType=='slow_alphaIon_betaIon') {
      unlink(paste0('survival-data-', particle), recursive = TRUE)
    }
  }
  
  return(out.df)
}


#' Evaluation of ion alpha and beta LQ parameters
#'
#' The evaluation is performed for a monoenergetic ion. It uses a C++ implementation for the evaluation. It accepts a sequence of parameters and ion specifications stored in an input data.frame.
#'
#' @param model the name of the model (options: 'MKM', 'LEMI', 'LEMII', 'LEMIII')
#' @param parameters input data.frame containing the cell name (column 'cell', optional) "biological" parameters associated to a specific biological tissue/model, and the specification of the primary particles (ion type and energy and/or LET). All the parameters have to be put in a dataframe. Each row of the dataframe represent the irradiation of a specific tissue/model with a specific particle. In the case of the MKM the biological parameters
#' \itemize{
#'   \item{alpha0} LQ alpha parameter of the reference radiation [Gy^-1]
#'   \item{beta0} LQ beta parameter of the reference radiation [Gy^-2]
#'   \item{rN} cell nucleus radius [um]
#'   \item{rd} domain radius [um]
#' }
#' In the case of LEMI, LEMII, LEMIII:
#' \itemize{
#'   \item{alpha0} LQ alpha parameter of the reference radiation [Gy^-1]
#'   \item{beta0} LQ beta parameter of the reference radiation [Gy^-2]
#'   \item{rN} cell nucleus radius [um]
#'   \item{Dt} threshold dose [Gy]
#' }
#' The ion is specified with the columns:
#' \itemize{
#'   \item{particle} ion type: 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F' or 'Ne'
#'   \item{energy} Energy of the particle in MeV/u
#'   \item{let} LET of the particle in keV/um. It is used only if the energy is not specified.
#' }
#' #' @family LEM/MKM Models
#' @export
#' @import data.table
alpha.beta.ion.sequence <- function(model='MKM',
                                    parameters=data.frame(cell='test', alpha0=0.1295, beta0=0.03085, rN=4, rd=0.31,
                                                          particle='H', energy=50, let=NA),
                                    calculusType='rapidMKM', precision=0.5,
                                    ignore.stdout=TRUE, ignore.stderr=TRUE,
                                    remove.temp.files=TRUE)
{
  
  #parameters <- data.table(parameters)
  parameters.column.index <- which(names(parameters) %in% c('alpha0', 'beta0', 'rN', 'rd', 'Dt'))
  if(length(parameters.column.index)!=4) {stop('alpha.beta.ion.sequence: error in the definition of the parameters')}
  
  # fa diventare data.table....
  #if('data.table' %in% class(parameters)) {
  #  message('casting in data.table')
  #  parameters <- as.data.table(parameters)
  #}
  
  # aggiunge identificativo unico per i parametri radiobiologici
  
  #parameters$id.par <- interaction(parameters[, parameters.column.index])
  
  # ho perso le chiavi
  #setkey(parameters, id.par)
  
  # dataframe di output
  alpha.beta.out <- NULL
  
  # loop su id.par
  Npar <- nrow(parameters)

  for(id in 1:Npar) {
    if(nrow(parameters)==0) {break}
    par <-  parameters[1, parameters.column.index]
    index.par <- parameters[,parameters.column.index[1]] == par[1,1] &
      parameters[,parameters.column.index[2]] == par[1,2] &
      parameters[,parameters.column.index[3]] == par[1,3] &
      parameters[,parameters.column.index[4]] == par[1,4]
    parameters.m <- parameters[index.par, ]
    
    cell <- parameters.m$cell[1]
    message('cell: ', cell)
    print(par)
    #print(summary(parameters.m))
    
    # loop su particles
    particles <- unique(parameters.m$particle)
    for(ip in 1:length(particles)) {
      message('particle: ', particles[ip])
      parameters.p <- subset(parameters.m, particle==particles[ip])
      
      energies <- parameters.p$energy
      lets <- parameters.p$let
      
      alpha.beta.tmp <- alpha.beta.ion.range(model=model, model.parameters=par, cell.name=cell, particle=particles[ip], energies=energies, lets=lets, calculusType=calculusType, precision=precision, ignore.stdout=ignore.stdout, ignore.stderr=ignore.stderr, remove.temp.files=remove.temp.files)
      alpha.beta.out <- rbind(alpha.beta.out, alpha.beta.tmp)
    }
    
    parameters <- parameters[-which(index.par),]
    
  }
  
  return(alpha.beta.out)  
}

#' @export
#' @import data.table
alpha.beta.ion.sequence2 <- function(model='MKM',
                                    parameters=data.frame(cell='test', alpha0=0.1295, beta0=0.03085, rN=4, rd=0.31,
                                                          particle='H', energy=50, let=NA),
                                    calculusType='rapidMKM', precision=0.5,
                                    ignore.stdout=TRUE, ignore.stderr=TRUE,
                                    remove.temp.files=TRUE)
{
  
  #parameters <- data.table(parameters)
  parameters.column.index <- which(names(parameters) %in% c('alpha0', 'beta0', 'rN', 'rd', 'Dt'))
  if(length(parameters.column.index)!=4) {stop('alpha.beta.ion.sequence: error in the definition of the parameters')}
  
  # fa diventare data.table....
  #if('data.table' %in% class(parameters)) {
  #  message('casting in data.table')
  #  parameters <- as.data.table(parameters)
  #}
  
  # aggiunge identificativo unico per i parametri radiobiologici

  #parameters$id.par <- interaction(parameters[, parameters.column.index])
  
  # ho perso le chiavi
  #setkey(parameters, id.par)
  
  # dataframe di output
  alpha.beta.out <- NULL
  
  # loop su id.par
  message('identifying models...')
  id.pars <- unique(parameters[, parameters.column.index])
  message(nrow(id.pars), ' different models')
  for(id in 1:nrow(id.pars)) {
    parameters.m <- parameters[parameters[,parameters.column.index[1]] == id.pars[id,1] &
                                  parameters[,parameters.column.index[2]] == id.pars[id,2] &
                                  parameters[,parameters.column.index[3]] == id.pars[id,3] &
                                  parameters[,parameters.column.index[4]] == id.pars[id,4], ]
    
    biological.parameters <- id.pars[id,]
    cell <- parameters.m$cell[1]
    message('cell: ', cell)
    print(biological.parameters)
    
    # loop su particles
    particles <- unique(parameters.m$particle)
    for(ip in 1:length(particles)) {
      message('particle: ', particles[ip])
      parameters.p <- subset(parameters.m, particle==particles[ip])
      
      energies <- parameters.p$energy
      lets <- parameters.p$let
      
      alpha.beta.tmp <- alpha.beta.ion.range(model=model, model.parameters=biological.parameters, cell.name=cell, particle=particles[ip], energies=energies, lets=lets, calculusType=calculusType, precision=precision, ignore.stdout=ignore.stdout, ignore.stderr=ignore.stderr, remove.temp.files=remove.temp.files)
      alpha.beta.out <- rbind(alpha.beta.out, alpha.beta.tmp)
    }
  }
  
  return(alpha.beta.out) 
}

#' Evaluate alpha function (MKM)
#'
#' Returns a function alpha(LET) or alpha(specific energy), using the MKM model. The functions is an interpolation functions that uses a set of alpha values evaluated over the defined set of LETs (or specific energies).
#'
#' @param alphaX,betaX,rN,rd The radiobiological parameter of the MKM.
#' @param cellType Optionally you can associate a cell name tag (I have to figure yet what use to do with that).
#' @param particleType The particle specification ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F')
#' @param energies,lets A vector of specific energies (MeV/u) or LETs (keV/um).
#' @param use.limits Filters the energies or the LETs values to rule out values too high or too low.
#'
#' @family LEM/MKM Models
#' @export
#'
alpha.fun.mkm <- function(alphaX=0.1295, betaX=0.03085, rN=4, rd=0.31,
                          cellType=NULL,
                          particleType='H',
                          energies=NULL, lets=NULL,
                          use.limits=FALSE, outmessages=FALSE)
{
  model <- 'MKM'
  calculusType <- 'rapidMKM'
  precision <- 1
  
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
  s.args <- paste(cellType, model, calculusType, a, b, n, d, particleType, precision, energyType, e)
  #cmd <- paste('.', lem.setenv, '; survival_alpha_beta_parameter_study', s.args)
  cmd <- paste('survival', s.args)
  #print(cmd)
  
  if(outmessages) {ignore.stdout=FALSE; ignore.stderr=FALSE} else {ignore.stdout=TRUE; ignore.stderr=TRUE}
  
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

