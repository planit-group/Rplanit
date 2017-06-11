# LUT --------------------------------------------------------------------------


#' Read LUT (4d)
#'
#' Read pencil-LUT.
#' @param lut.name the name (prefix of the file name) of the lut (i.e. not including the file name extension.)
#' @param dataframe return a dataframe.
#' @family LUT
#' @export
read.lut <- function(lut.name, dataframe=FALSE)
{

  lut.file <- paste(lut.name, '.lut.4d', sep='')
  #lut.file <- lut.name

  # legge lutMeta. Serve per identificare i tessuti.
  N.tissues <- read.table(file = paste0(lut.name, '.lutMeta'), skip = 4, nrows = 1)
  tissues <- as.vector(read.table(file = paste0(lut.name, '.lutMeta'), skip = 5, nrows = as.numeric(N.tissues))[1])
  message('Found ', N.tissues, ' tissues...')
  print(tissues)

  # legge file 4d
  message('reading lut:', lut.file)

  # apri connsessione
  file.4d <- file(lut.file, "rb") # read binary

  # leggi l'header
  message('reading header...')
  myline <- readLines(file.4d, n=10) # legge le prime 10 linee, ogni linea e' un elemento del vettore

  # parsing dell'header
  # splitta la stringa in sottostringhe delimitate da uno o piu' spazi (regexpr: " +")

  # NominalEnergy[MeV/u]
  myline.splitted <- unlist(strsplit(myline[1], ' +'))
  NE <- as.numeric(myline.splitted[1])
  E.name <- myline.splitted[2]
  E.type <- myline.splitted[3]
  cat(NE, E.name, E.type, '\n')
  myline.splitted <- unlist(strsplit(myline[2], ' +'))
  E <- as.numeric(myline.splitted)

  # X[mm]
  myline.splitted <- unlist(strsplit(myline[3], ' +'))
  NX <- as.numeric(myline.splitted[1])
  X.name <- myline.splitted[2]
  X.type <- myline.splitted[3]
  cat(NX, X.name, X.type, '\n')
  myline.splitted <- unlist(strsplit(myline[4], ' +'))
  X <- as.numeric(myline.splitted)

  # Y[mm]
  myline.splitted <- unlist(strsplit(myline[5], ' +'))
  NY <- as.numeric(myline.splitted[1])
  Y.name <- myline.splitted[2]
  Y.type <- myline.splitted[3]
  cat(NY, Y.name, Y.type, '\n')
  myline.splitted <- unlist(strsplit(myline[6], ' +'))
  Y <- as.numeric(myline.splitted)

  # Z[mm]
  myline.splitted <- unlist(strsplit(myline[7], ' +'))
  NZ <- as.numeric(myline.splitted[1])
  Z.name <- myline.splitted[2]
  Z.type <- myline.splitted[3]
  cat(NZ, Z.name, Z.type, '\n')
  myline.splitted <- unlist(strsplit(myline[8], ' +'))
  Z <- as.numeric(myline.splitted)

  # Variables
  myline.splitted <- unlist(strsplit(myline[9], ' +'))
  NV <- as.numeric(myline.splitted[1])
  variables <- myline.splitted[2:length(myline.splitted)]

  # Format
  myline.splitted <- unlist(strsplit(myline[10], ' +'))
  values.format <- myline.splitted[1]

  Ntot <- NE * NX * NY * NZ * NV
  cat('number of voxels:', NV, 'x', NE, 'x', NX, 'x', NY, 'x', NZ, '=', Ntot, '\n')
  cat('dependent variables:', variables, '\n')
  cat('values format:', values.format,'\n')

  # trasforma intervalli in coordinate puntuali
  if(E.type=='INTERVAL') {E <- (E[1:NE] + E[2:(NE+1)])/2}
  if(X.type=='INTERVAL') {X <- (X[1:NX] + X[2:(NX+1)])/2}
  if(Y.type=='INTERVAL') {Y <- (Y[1:NY] + Y[2:(NY+1)])/2}
  if(Z.type=='INTERVAL') {Z <- (Z[1:NZ] + Z[2:(NZ+1)])/2}

  # crea e legge array (4d+1)
  if(values.format=='BINARY') {
    message('reading binary data...')
    Values.4d <- array(readBin(file.4d, numeric(), Ntot), dim=c(NV, NE, NX, NY, NZ))
  } else {
    message('reading ascii data...')
    myvector <- as.numeric(unlist(strsplit(readLines(file.4d), ' +')))
    #print(myvector)
    Values.4d <- array(myvector, dim=c(NV, NE, NX, NY, NZ))
    rm(myvector)
  }

  # chiudi connessione
  close(file.4d)

  #return(Values.4d)

  # z.BP per calcolare Z assoluti
  z.BP <- get.zBP.fun(lut.name)
  #z.BP <- 1

  # aggiunge tag per i diversi tessuti
  # alpha
  v.tissues <- which(variables=='LethAlphaPerEvent[1/primary]')
  if(length(v.tissues)>0) {
    for(i in seq_along(v.tissues)) {
      variables[v.tissues[i]] <- paste(variables[v.tissues[i]], tissues[i,1], sep=' ')
    }
  }
  # beta
  v.tissues <- which(variables=='SqrtLethBetaPerEvent[1/primary]')
  if(length(v.tissues)>0) {
    for(i in seq_along(v.tissues)) {
      variables[v.tissues[i]] <- paste(variables[v.tissues[i]], tissues[i,1], sep=' ')
    }
  }
  print(variables)



  # ritorna lista+array
  if(!dataframe) {
    return(list(values=Values.4d, E=E, x=X, y=Y, z=Z, z.BP=z.BP, NE=NE, Nx=NX, Ny=NY, Nz=NZ, Nv=NV, variables=variables))
  }

  # crea dataframe strutturato
  message('melting data into data.frame...')
  Values <- melt(Values.4d)
  rm(Values.4d)
  if (NV==1) {
    Values <- cbind(1, Values)
  }
  names(Values) <- c('variable', 'E', 'x', 'y', 'z', 'value')

  Values$variable <- as.factor(sapply(Values$variable, function(v) variables[v]))

  Values$E <- (sapply(Values$E, function(v) E[v]))
  Values$x <- (sapply(Values$x, function(v) X[v]))
  Values$y <- (sapply(Values$y, function(v) Y[v]))
  Values$z <- (sapply(Values$z, function(v) Z[v]))
  Values$z.Bp <- z.BP(Values$E)
  Values$za <- Values$z * Values$z.Bp

  return(Values)
}

#' Write LUT meta file
#' @param lut.name The name (file prefix) of the LUT
#' @param physical.variables Names of the physical variables
#' @param biological.variables Names of the biological variables
#' @param model.names Names of the biological models
#' @param alphaX,betaX vectors of the LQ parameters for the reference radiation (used to compute the RBE.)
#' @family LUT
#' @export
write.lut.meta <- function(lut.name,
                           physical.variables='DosePerEvent[Gy/primary]',
                           biological.variables=c('LethAlphaPerEvent[1/primary]', 'SqrtLethBetaPerEvent[1/primary]'),
                           tissue.names='NULL', alphaX=1, betaX=1) {
  file.name <- paste0(lut.name, '.lutMeta')
  con <- file(file.name, "w") # open for writing in text mode
  
  # optics
  writeLines(paste0(LUT.name, ".optics\n"), con=con, sep='') # Name of beam-line-optics file
  writeLines("0.0   0.0   0.0\n", con=con, sep='') # Deflection along X, deflection along Y, phantom entrance along Z
  
  # physical variables
  physical.text <- paste(length(physical.variables), paste0(physical.variables, collapse = ' '), '\n')
  writeLines(physical.text, con=con, sep='') # number of physical values, list of physical values tags
  
  # biological variables
  if(length(biological.variables)>0) {
    biological.text <- paste(length(biological.variables), paste0(biological.variables, collapse = ' '), '\n')
    writeLines(biological.text, con=con, sep='') # number of biological values, list of biological values
    writeLines(paste0(length(tissue.names), '\n'), con=con, sep='')
    for(im in seq_along(tissue.names)) {
      writeLines(paste0(tissue.names[im], ' ', alphaX[im], ' ', betaX[im], '\n'), con=con, sep='')
    }
  }
  
  # files
  writeLines(paste0(LUT.name, ".zBraggPeaks.1d\n"), con=con, sep='') # Name of Bragg-peak-depths file
  writeLines(paste0(LUT.name, ".lut.4d\n"), con=con, sep='') # Name of LUT file
  
  close(con)
}

#' Write LUT optics (meta) file
#' @param lut.name The name (file prefix) of the LUT
#' @param x.magnet.distance,y.magnet.distance Distances of the scanning magnets.
#' @param particle The primary particle.
#' @family LUT
#' @export
write.lut.optics <- function(lut.name, x.magnet.distance=6000, y.magnet.distance=6000, particle='1H') {
  file.name <- paste0(lut.name, '.optics')
  con <- file(file.name, "w") # open for writing in text mode
  writeLines(particle, con=con) # particle type
  writeLines(paste(x.magnet.distance, y.magnet.distance), con=con) # virtual source vsad_x vsad_y
  writeLines(paste0(LUT.name, ".sub.lutMeta"), con=con) # Name of sub-beams meta file
  writeLines(paste0(LUT.name, ".xPhaseSpace.3d"), con=con) # Name of x-phase-space file
  writeLines(paste0(LUT.name, ".yPhaseSpace.3d"), con=con) # Name of y-phase-space file
  writeLines(paste0(LUT.name, ".energyPhaseSpace.1d"), con=con) # Name of energy-phase-space file
  writeLines(paste0(LUT.name, ".shiftsPhaseSpace.1d"), con=con) # Name of z-shifts-phase-space file
  close(con)
}

#' Write LUT Bragg peaks depths file
#' @param lut.name The name (file prefix) of the LUT
#' @param energies Vector of energies.
#' @param zBP Corresponding vector of Bragg peack depths [mm].
#' @family LUT
#' @export
write.lut.zBP <- function(lut.name, energies, zBP) {
  file.name <- paste0(lut.name, '.zBraggPeaks.1d')
  con <- file(file.name, "w") # open for writing in text mode
  writeLines(paste0(length(energies), ' Energy[MeV/u] POINTWISE'), con=con) # Number of energies
  writeLines(paste(energies, collapse = ' '), con=con) # List of energies
  writeLines('1 BraggPeakDepths[mm]', con=con) # number of values (1: bragg peak depths)
  writeLines('ASCII', con=con) # formato
  writeLines(as.character(zBP), con=con)
  close(con)
}

#' Write LUT 4d data array
#' @param lut.name The name (file prefix) of the LUT
#' @param lut the LUT object containing the 4d data array
#' @family LUT
#' @export
write.lut.4d <- function(lut.name, lut) {
  # la lut è sostanzialmente un oggetto values 4d:
  # le dimensioni dell'oggetto 4d sono: NV, NE, NX, NY, NZ
  # Le coordinate Z sono normalizzate e alla lut è agganciata anche una funzione z.BP(z), per ottenere i valori assoluti di z.
  
  # nota: nell'oggetto lut alphaD e sqrt(beta)D sono differenziati con '*_n', mentre nel
  # file da salvare no...
  
  variables <- gsub('_[0-9]|_[1-9][0-9]', '', lut$variables)
  # print(variables)
  
  
  # 4d file
  
  # scrive header
  file.name <- paste0(lut.name, '.lut.4d')
  con <- file(file.name, "w") # open for writing in text mode
  
  writeLines(paste(lut$NE, "NominalEnergy[MeV/u] POINTWISE"), con=con) # Number of energies
  writeLines(paste(lut$E, collapse=' '), con=con)
  writeLines(paste(lut$Nx, "X[mm] POINTWISE"), con=con)
  writeLines(paste(lut$x, collapse=' '), con=con)
  writeLines(paste(lut$Ny, "Y[mm] POINTWISE"), con=con)
  writeLines(paste(lut$y, collapse=' '), con=con)
  writeLines(paste(lut$Nz, "NormalizedZ POINTWISE"), con=con)
  writeLines(paste(lut$z, collapse=' '), con=con)
  writeLines(paste(lut$Nv, paste(variables, collapse=' ')), con=con) # Nv = Number of physical quantities + (number of biological quantities * number of tissues)
  writeLines('BINARY', con=con)
  close(con)
  
  # scrive dati binari
  con <- file(file.name, "ab") # open for writing in text mode
  writeBin(as.double(lut[['values']]), con=con)
  close(con)
}

#' Write LUT phase-space file
#' 
#' Note: this implementations assumes a "complete" pencilbeam approach, without sub-pencilbeams. Hence the phase space is essentially
#' a place holder.
#' @param lut.name The name (file prefix) of the LUT
#' @param lut the LUT object containing the 4d data array
#' @family LUT
#' @export
write.lut.phasespace <- function(lut.name,
                                 nominalEnergy=NULL,
                                 meanEnergy=NULL, sigmaEnergy=0,
                                 sigmaX=4, sigmaY=4, sigmaX1=0, sigmaY1=0, rhoXX1=0, rhoYY1=0) {
  if(is.null(nominalEnergy)) {nominalEnergy <- c(1,1e3)}
  NE <- length(nominalEnergy)
  
  if(is.null(meanEnergy)) {meanEnergy <- nominalEnergy}
  
  if(length(sigmaX)==1) {sigmaX <- rep(sigmaX, NE)}
  if(length(sigmaY)==1) {sigmaY <- rep(sigmaY, NE)}
  if(length(sigmaX1)==1) {sigmaX1 <- rep(sigmaX1, NE)}
  if(length(sigmaY1)==1) {sigmaY1 <- rep(sigmaY1, NE)}
  if(length(rhoXX1)==1) {rhoXX1 <- rep(rhoXX1, NE)}
  if(length(rhoYY1)==1) {rhoYY1 <- rep(rhoYY1, NE)}
  if(length(sigmaEnergy)==1) {sigmaEnergy <- rep(sigmaEnergy, NE)}
  
  # nota: i valori sigmaX,Y e sigmaX1,Y1 dipendono in generale dall'energia e dalla deflessione.
  # Qui per semplicfictà si assume che non ci sia nessuna dipendenza in funzione della defessione
  # Sembra che invece lo spazio delle fasi dell'energia non dipenda dalla deflessione... boh?
  
  # x phasespace
  file.name <- paste0(lut.name, '.xPhaseSpace.3d')
  con <- file(file.name, "w")
  writeLines(paste(NE, 'NominalEnergy[MeV/u]  POINTWISE'), con=con)
  writeLines(paste(nominalEnergy, collapse=' '), con=con)
  writeLines('2 DeflectionX[mm] POINTWISE', con=con)
  writeLines('-1e9 1e9', con=con)
  writeLines('2 DeflectionY[mm] POINTWISE', con=con)
  writeLines('-1e9 1e9', con=con)
  writeLines("3 SigmaX[mm] SigmaX' RhoXX'", con=con)
  writeLines('ASCII', con=con)
  for(ii in 1:4) { # ripete valori identificativi al variare delle deflessioni
    writeLines(paste(sigmaX, sigmaX1, rhoXX1), con=con)
  }
  close(con)
  
  # y phasespace
  file.name <- paste0(lut.name, '.yPhaseSpace.3d')
  con <- file(file.name, "w")
  writeLines(paste(NE, 'NominalEnergy[MeV/u]  POINTWISE'), con=con)
  writeLines(paste(nominalEnergy, collapse=' '), con=con)
  writeLines('2 DeflectionX[mm] POINTWISE', con=con)
  writeLines('-1e9 1e9', con=con)
  writeLines('2 DeflectionY[mm] POINTWISE', con=con)
  writeLines('-1e9 1e9', con=con)
  writeLines("3 SigmaY[mm] SigmaY' RhoYY'", con=con)
  writeLines('ASCII', con=con)
  for(ii in 1:4) { # ripete valori identivi al variare delle deflessioni
    writeLines(paste(sigmaY, sigmaY1, rhoYY1), con=con)
  }
  close(con)
  
  # energy phasespace
  file.name <- paste0(lut.name, '.energyPhaseSpace.1d')
  con <- file(file.name, "w")
  writeLines(paste(NE, 'NominalEnergy[MeV/u]  POINTWISE'), con=con)
  writeLines(paste(nominalEnergy, collapse=' '), con=con)
  writeLines('2 MeanEnergy[MeV/u] SigmaEnergy[MeV/u]', con=con)
  writeLines('ASCII', con=con)
  writeLines(paste(meanEnergy, sigmaEnergy), con=con)
  close(con)
  
  # shift phasespace
  file.name <- paste0(lut.name, '.shiftsPhaseSpace.1d')
  con <- file(file.name, "w")
  writeLines(paste(NE, 'NominalEnergy[MeV/u]  POINTWISE'), con=con)
  writeLines(paste(nominalEnergy, collapse=' '), con=con)
  writeLines('2 ShiftZ_1[mm] Weight_1', con=con) # ogni shift è una variabile a parte... ma che cazzo...
  writeLines('ASCII', con=con)
  writeLines(paste(c(0,0), c(1,1)), con=con)
  close(con)
}

#' Return a function z.BP(E)
#'
#' z.BP is the depth of the Bragg peak as a function of Energy (E). Normally the energy is intended as the specific energy (MeV/u). The function is obtained by interpolating the data contained in the file <lut.name>.zBraggPeaks.1d.
#' @param lut.name The name of the lut.
#' @param EBP Return E(z.BP) instead of z.BP(E) if true.
#' @family LUT
#' @export
get.zBP.fun <- function(lut.name, EBP=FALSE)
{

  zBP.file <- paste(lut.name, '.zBraggPeaks.1d', sep='')

  # legge file 1d
  message('reading lut:', zBP.file)
  file.1d <- file(zBP.file, "rb")
  myline <- readLines(file.1d, n=4)
  # Energy[MeV/u]
  myline.splitted <- unlist(strsplit(myline[1], ' +'))
  NE.BP <- as.numeric(myline.splitted[1])
  E.BP.name <- myline.splitted[2]
  E.BP.type <- myline.splitted[3]
  cat(NE.BP, E.BP.name, E.BP.type, '\n')
  myline.splitted <- unlist(strsplit(myline[2], ' +'))
  E.BP <- as.numeric(myline.splitted)

  # Format
  myline.splitted <- unlist(strsplit(myline[4], ' +'))
  values.format <- myline.splitted[1]

  if(values.format=='BINARY') {
    message('reading binary data...')
    z.BP <- array(readBin(file.1d, numeric(), Ntot), dim=c(NE.BP))
  } else {
    message('reading ascii data...')
    z.BP <- as.numeric(unlist(strsplit(readLines(file.1d), ' +')))
  }

  close(file.1d)

   if(EBP) {
     return(approxfun(z.BP, E.BP))
   } else {
     return(approxfun(E.BP, z.BP))
   }
}


#' Add alpha values to lUT
#'
#' Add variable Alpha[Gy^(-1)] = LethAlphaPerEvent / DosePerEvent
#'
#' @param lut the lut object
#' @param tissue.number the number of the tissue stored in the LUT (start from 1)
#' @family LUT
#' @export
lut.add.alpha <- function(lut, tissue.number=1, tissue=NULL)
{
  if(!is.null(tissue)) {
    variable.alpha <- paste('LethAlphaPerEvent[1/primary]', tissue, sep=' ')
  } else {
    variable.alpha <- lut$variables[which(lut$variables=='LethAlphaPerEvent[1/primary]')[tissue.number]]
    tissue <- tissue.number
  }
  message('evaluating alpha... ')
  if(class(lut)=='data.frame'){
    a <- subset(lut, variable==variable.alpha)
    d <- subset(lut, variable=='DosePerEvent[Gy/primary]')
    a$value <- a$value / d$value
    a$variable <- 'Alpha[Gy^(-1)]'
    return(rbind(lut, a))
  } else {
    Na <- which(lut$variables==variable.alpha)
    Nd <- which(lut$variables=='DosePerEvent[Gy/primary]')
    a <- lut$values[Na,,,,, drop=FALSE]
    d <- lut$values[Nd,,,,, drop=FALSE]
    a <- a / d
    #plot(a)
    av <- aperm(lut$values, c(2,3,4,5,1))
    a <- aperm(a, c(2,3,4,5,1))
    av <- array(c(av,a), c(lut$NE, lut$Nx, lut$Ny, lut$Nz, lut$Nv+1) )
    lut$values <- aperm(av, c(5,1,2,3,4))
    lut$variables[lut$Nv+1] <- paste('Alpha[Gy^(-1)]', tissue)
    lut$Nv <- lut$Nv+1
    return(lut)
  }
}


#' aggiungi beta alla lut
#'
#' la formula: LethAlphaPerEvent / DosePerEvent
#' @family LUT
#' @export
lut.add.beta <- function(lut)
{
  message('evaluating beta... ')
  if(class(lut)=='data.frame'){
    b <- subset(lut, variable=='SqrtLethBetaPerEvent[1/primary]')
    d <- subset(lut, variable=='DosePerEvent[Gy/primary]')
    b$value <- b$value / d$value
    b$value <- b$value * b$value
    b$variable <- 'Beta[Gy^(-2)]'
    return(rbind(lut, b))
  } else {
    message('lut.add.beta for array+list not yet implemented...')
    return()
  }
}


#' Aggiunge RBE.alpha alla LUT
#' @family LUT
#' @export
lut.add.rbe.alpha <- function(lut, alphaX)
{
  message('evaluating RBE.alpha... ')
  if(class(lut)=='data.frame'){
    a <- subset(lut, variable=='LethAlphaPerEvent[1/primary]')
    d <- subset(lut, variable=='DosePerEvent[Gy/primary]')
    a$value <- a$value / d$value / alphaX
    a$variable <- 'RBE.alpha'
    return(rbind(lut, a))
  } else {
    Na <- which(lut$variables=='LethAlphaPerEvent[1/primary]')
    Nd <- which(lut$variables=='DosePerEvent[Gy/primary]')
    a <- lut$values[Na,,,,, drop=FALSE]
    d <- lut$values[Nd,,,,, drop=FALSE]
    a <- a / d / alphaX
    av <- aperm(lut$values, c(2,3,4,5,1))
    a <- aperm(a, c(2,3,4,5,1))
    av <- array(c(av,a), c(lut$NE, lut$Nx, lut$Ny, lut$Nz, lut$Nv+1) )
    lut$values <- aperm(av, c(5,1,2,3,4))
    lut$variables[lut$Nv+1] <- 'RBE.alpha'
    lut$Nv <- lut$Nv+1
    return(lut)
  }
}

#' integra la lut radiale
#' @family LUT
#' @export
lut.radial.integral <- function(lut, r.cut=Inf)
{
  message('evaluating radial lut integral...')
  if(class(lut)!='data.frame') {
    message('the lut is not a data.frame...')
    return()
  }
  if(r.cut<Inf) {lut <- subset(lut, x<r.cut)}
  r <- sort(unique(lut$x))
  r.a <- r*r
  r.w <- c(0, diff(r.a))
  r.w.fun <- approxfun(r, r.w)
  lut$r.w <- r.w.fun(lut$x)
  lut.int <- aggregate(data.frame(value=lut$value*lut$r.w),
                       list(E=lut$E, z=lut$z, za=lut$za, z.Bp=lut$z.Bp, variable=lut$variable),
                       mean)
  return(lut.int)
}


#' salva la lut a partire dai dati in memoria
#' @family LUT
#' @export
write.lut <- function(lut.array, variables=NULL, E=NULL, x=NULL, y=NULL, zn=NULL, filename='lut.4d', binary=FALSE)
{

  Ne <- length(E)
  Nx <- length(x)
  Ny <- length(y)
  Nzn <- length(zn)
  Nv <- length(variables)

  my.file <- file(filename, "wt")

  writeLines(paste(Ne, ' NominalEnergy[MeV/u]   POINTWISE\n'), con=my.file, sep='')
  writeLines(paste(E, con=my.file, collapse=' ')); writeLines('\n', con=my.file)

  writeLines(paste(Nx, ' X[mm]   POINTWISE\n'), con=my.file, sep='')
  writeLines(paste(x, con=my.file, collapse=' ')); writeLines('\n', con=my.file)

  writeLines(paste(Ny, ' Y[mm]   POINTWISE\n'), con=my.file, sep='')
  writeLines(paste(y, con=my.file, collapse=' ')); writeLines('\n', con=my.file)

  writeLines(paste(Nzn, ' NormalizedZ   POINTWISE\n'), con=my.file, sep='')
  writeLines(paste(zn, con=my.file, collapse=' ')); writeLines('\n', con=my.file)

  writeLines(paste(Nv, ' '), con=my.file)
  writeLines(paste(variables, con=my.file, collapse=' ')); writeLines('\n', con=my.file)

  close(my.file)

  if(binary) {
    my.file <- file(filename, "ab")
    writeBin(lut.array, my.file)
  } else {
    message('ASCII not implemented.')
  }

  return()
}


# BEAM_LUT ---------------------------------------------------------------------

#' Get beamLUTs
#'
#' Get the evaluated beamLUTs from the plan
#' @param plan The plan object.
#' @param concatenate use alternative faster reading approach (parameter passed to read.beamLUT(), may not work on windows and osx.)
#' @export
#' @family BeamLUT
get.beamLUTs <- function(plan, concatenate=FALSE) {
  if(is.null(plan[['outputBeamLUTFile']])) {
    stop('beamLUTs not present in plan ', plan[['name']], '. To evaluate the beamLUT set the plan with saveBeamLUTs=TRUE.')
  } else {
    values <- get.values(plan)
    return(read.beamLUTs(beamLUT.name = plan[['outputBeamLUTFile']], Nx = values$Nx, Ny = values$Ny, Nz = values$Nz, concatenate = concatenate))
  }
}

#' Evaluate BeamLUT (obsolete)
#'
#' Evaluate from scratch the beamLUT of the selected beam(s) for the specific plan.
#' The beam(s) can be one or more (by using a vector of indices). If the beams are more than one, the corresponding beamLUT is given by the net contribution of all the specified beams
#' Optionally it is possible to specify also the fluences (i.e. the number of particles) for each beam.
#' @param beam.index index of the beam (it can be a vector if indices)
#' @param fluence the fluence (number of particles) for the beam(s). It can be a vector
#' @param plan the plan object
#' @param variable string specifing the variable to be evaluated. It can be a vecor of variables
#' @param temp.name temporary prefix of the plan name used during the evaluation
#' @param remove.temp TRUE to remove temporary files after the evaluation
#' @return a values object
#' @family BeamLUT
#' @export
get.values.for.beam <- function(beams=NULL, beam.index, fluence=NULL, plan, variables='Dose[Gy]', temp.name='temp', remove.temp=TRUE)
{
  # recupera e seleziona beams
  if(is.null(beams)) {
    beams <- get.beams(plan)
    beams <- beams[beam.index, ]
  }

  temp.beamfile <- paste(temp.name, '.beams', sep='')
  if(!is.null(fluence)) {beams$fluence=fluence; print(beams)} # sovrascrive fluenze
  write.beams(beams, temp.name)

  # crea piano pre la selezione dei beam e calcola forward.planning
  plan.b <- plan
  plan.b$name <- temp.name
  plan.b['beams'] <- list(NULL) # disinnesca i beam (se per caso ci sono...)
  plan.b$inputBeamsFile <- temp.beamfile
  plan.b$computingValues <- paste(variables, collapse=' ')
  plan.b.out <- run.dek.forward(plan.b, outmessages=FALSE)

  # recupera values
  vb <- get.values(plan.b.out)

  # rimuovi file temporanei
  if(remove.temp) {
    unlink(temp.beamfile)
    unlink(temp.name, recursive=TRUE)
  }

  # ritorna values
  return(vb)
}


#' scrive le matrici V_{b,v} (beamLUT) (obsolete)
#'
#' nota: è usata la "threshold.variable" (normalmente la dose) per fissare la soglia
#' le lut sono normalizzate per singolo primario
#' le lut sono calcolate su tutto il volume usato in "plan"
#'
#' @family BeamLUT
#' @export
write.beamLUT <- function(plan, threshold=0, threshold.variable='Dose[Gy]', variable='Dose[Gy]', file.name)
{

  beams <- get.beams(plan)
  beams$beamID <- 1:nrow(beams)
  #beams <- subset(beams, fluence>0)

  if(threshold.variable==variable) {variables=variable} else {variables=c(threshold.variable, variable)}

  for(b in 1:nrow(beams)) {
    message('writing beam ', b)
    vb <- get.values.for.beam(beam.index=beams$beamID[b], fluence=1, plan=plan, variables=variables)
    df.s <- sparse.array.from.values(values=vb, variable=variable, threshold=threshold)
    #df.s$value <- df.s$value/beams$fluence[b]
    df.s$beamID <- beams$beamID[b]

      if(b==1) {append=FALSE} else {append=TRUE}
      write.table(df.s, file=file.name, col.names=FALSE, row.names=FALSE, append=append)

  }
}


#' Read a beamLUT (pure-dek)
#'
#' @param beamLUT.file the name of the beamLUT file.
#' @param Nx,Ny,Nz optional dimensions of the target 3D array (to evaluate an absolute voxelID)
#' @param dose.threshold additional dose threshold to filter elements with dose <= dose.threshold.
#' @return a data.table of the beamLUT
#' @family BeamLUT
#' @export
#' @importFrom data.table fread
read.beamLUT <- function(beamLUT.file, Nx=NULL, Ny=NULL, Nz=NULL, dose.threshold=0)
{
  beamLUT <- fread(beamLUT.file, skip = 1, header = TRUE, check.names=FALSE)
  
  # filtro dose.thresold
  if(dose.threshold > 0) {
    beamLUT <- subset(beamLUT, `DosePerEvent[Gy/primary]` > dose.threshold)
  }
  
  # beamID start from 1...
  beamLUT$beamID <- beamLUT$beamID + 1
  
  # calcola i voxelID...
  if(!is.null(Nx) & !is.null(Ny) & !is.null(Nz)) {
    beamLUT$voxelID <- beamLUT$i + beamLUT$j*Nx + beamLUT$k*(Nx*Ny) + 1
  }
  return(beamLUT)
}

#' Read the beamLUTs (pure-dek)
#'
#' @param beamLUT.name the prefix of the beamLUT files. The complete filenames are expected to be: <beamLUT.name>_<beamIndex>.beamLUT
#' @param Nx,Ny,Nz the dimensions of the target 3D array (to evaluate an absolute voxelID)
#' @param dose.threshold additional dose threshold to filter elements with dose <= dose.threshold.
#' @param concatenate uses alternative faster approach (may not work on windows and osx.)
#' @family BeamLUT
#' @export
#' @importFrom plyr ldply
#' @importFrom data.table fread setDT
read.beamLUTs <- function(beamLUT.name, Nx, Ny, Nz, dose.threshold=0, concatenate=FALSE)
{
  beamLUT.files <- Sys.glob(paste0(beamLUT.name, '*.beamLUT'))
  if(length(beamLUT.files)==0) {stop(paste0(beamLUT.name, '*.beamLUT not found!'))}
  Nb <- length(beamLUT.files)
  
  # ciclo preliminare per contare tutti i voxels
  total.voxels <- 0
  for(b in 1:Nb) {
    total.voxels <- total.voxels + read.table(file = beamLUT.files[b], header = FALSE, nrow = 1)[1]
  }
  message('found ', total.voxels, ' elements in beamLUTs (not filtered.)')
  
  if(concatenate) {
    # sistema con concatenazione...
    message('reading ', Nb, ' concatenated beamLUTs...')
    concantenate.cmd <- paste0('tail -q -n +3 ', beamLUT.name, '*.beamLUT')
    beamLUT <- fread(concantenate.cmd, header = FALSE)
    beamLUT.tmp <- read.beamLUT(beamLUT.files[1])
    names(beamLUT) <- names(beamLUT.tmp)
  } else {
    message('reading ', Nb, ' beamLUTs... ')
    beamLUT <- setDT(ldply(beamLUT.files, read.beamLUT, .progress='text'))
  }
  
  # filtro dose.thresold
  if(dose.threshold > 0) {
    beamLUT <- subset(beamLUT, `DosePerEvent[Gy/primary]` > dose.threshold)
  }
  
  # beamID start from 1...
  #beamLUT$beamID <- beamLUT$beamID + 1
  
  # calcola i voxelID... (nota: non viene calcolato in read.beamLUT)
  beamLUT$voxelID <- beamLUT$i + beamLUT$j*Nx + beamLUT$k*(Nx*Ny) + 1
  return(beamLUT)
}

#' Get Dose from beamLUTs (pure-dek)
#'
#' Get the net Dose distribution (values object or sparse array) from the beamLUTs. If coordinates x, y and z are specified it returns a complete values object
#' @param beamLUTs the beamLUTs data frame.
#' @param beams the beams data frame.
#' @param x,y,z the coordinates of the 3D array.
#' @param variable to override the default name of the variable ('Dose[Gy]').
#' @return a values object (if  x, y and z are specified), otherwise a sparse array.
#' @family BeamLUT
#' @export
get.dose.beamLUT <- function(beamLUTs, beams, x=NULL, y=NULL, z=NULL, variable='Dose[Gy]')
{

  # calcolo con aggregate
  dose.bl <- beamLUTs$`DosePerEvent[Gy/primary]` * beams$fluence[beamLUTs$beamID]
  dose.bl <- aggregate(list(dose=dose.bl), by = list(voxelID=beamLUTs$voxelID), sum)

  if(is.null(x) | is.null(y) | is.null(z)) { # matrice sparsa
    names(dose.bl) <- c('voxelID', variable)
    return(dose.bl)
  } else { # oggetto values completo
    dose <- array(0, dim=c(length(x), length(y), length(z)))
    dose[dose.bl$voxelID] = dose.bl$dose
    return(create.values(array.values = dose, variables = variable, x = x, y = y, z = z))
  }
}

#' Get dose averaged LET from beamLUTs (pure-dek)
#'
#' Get the net dose averged LET distribution (values object or sparse array) from the beamLUTs.
#' @param beamLUTs the beamLUTs data frame.
#' @param beams the beams data frame.
#' @param x,y,z the coordinates of the 3D array.
#' @param variable to override the default name of the variable ('DoseAveragedLET[keV/um]').
#' @return A values object (if  x, y and z are specified), otherwise a sparse array.
#' @family BeamLUT
#' @export
get.letd.beamLUT <- function(beamLUTs, beams, x=NULL, y=NULL, z=NULL, variable='DoseAveragedLET[keV/um]')
{
  endens.bl <- beamLUTs$`EnergyDensityPerEvent[keV/(um^3*primary)]` * beams$fluence[beamLUTs$beamID]
  letd.bl <- beamLUTs$`DoseWeightedEnergyDensityPerEvent[keV^2/(um^4*primary)]` * beams$fluence[beamLUTs$beamID]
  letd.bl <- aggregate(list(endens=endens.bl, letd=letd.bl), by = list(voxelID=beamLUTs$voxelID), sum)

  if(is.null(x) | is.null(y) | is.null(z)) {
    letd.bl <- data.frame(letd.bl$voxelID, letd.bl$letd/letd.bl$endens)
    names(letd.bl) <- c('voxelID', variable)
    return(letd.bl)
  } else {
    letd <- array(NA, dim=c(length(x), length(y), length(z)))
    letd[letd.bl$voxelID] = letd.bl$letd/letd.bl$endens
    return(create.values(array.values = letd, variables = variable, x = x, y = y, z = z))
  }
}

#' Get RBE from beamLUTs (pure-dek)
#'
#' Get the net RBE distribution (values object or sparse array) from the beamLUTs.
#' @param beamLUTs the beamLUTs data frame.
#' @param beams the beams data frame.
#' @param x,y,z the coordinates of the 3D array.
#' @param alphaX,betaX the Linear-Quadratic parameters of the reference radiation (they can be two arrays/vectors)
#' @param variable to override the default name of the variable ('RBE').
#' @return A values object (if  x, y and z are specified), otherwise a sparse array.
#' @family BeamLUT
#' @export
get.rbe.beamLUT <- function(beamLUTs, beams, x=NULL, y=NULL, z=NULL, variable='RBE', alphaX, betaX)
{


  dose.bl <- beamLUTs$`DosePerEvent[Gy/primary]` * beams$fluence[beamLUTs$beamID]
  alpha.bl <- beamLUTs$`LethAlphaPerEvent[1/primary]` * beams$fluence[beamLUTs$beamID]
  beta.bl <- beamLUTs$`SqrtLethBetaPerEvent[1/primary]` * beams$fluence[beamLUTs$beamID]

  rbe.bl <- aggregate(list(dose=dose.bl, alpha=alpha.bl, beta=beta.bl), by = list(voxelID=beamLUTs$voxelID), sum)
  i.dose <- rbe.bl$dose>0
  rbe.bl$alpha <- rbe.bl$alpha/rbe.bl$dose
  rbe.bl$beta <- rbe.bl$beta^2/rbe.bl$dose^2

  if(is.null(x) | is.null(y) | is.null(z)) {
    rbe.bl <- data.frame(rbe.bl$voxelID, rbe.evaluate(alpha=rbe.bl$alpha, beta=rbe.bl$beta, dose=rbe.bl$dose, alphaX=alphaX, betaX=betaX))
    names(rbe.bl) <- c('voxelID', variable)
    return(rbe.bl)
  } else {
    rbe <- array(NA, dim=c(length(x), length(y), length(z)))
    rbe[rbe.bl$voxelID] = rbe.evaluate(alpha=rbe.bl$alpha, beta=rbe.bl$beta, dose=rbe.bl$dose, alphaX=alphaX, betaX=betaX)
    return(create.values(array.values = rbe, variables = variable, x = x, y = y, z = z))
  }
}

#' Get alpha/beta from beamLUTs (pure-dek)
#'
#' Get the net alpha and beta Linear-Quadratic parameters distribution (values object or sparse array) from the beamLUTs.
#' @param beamLUTs the beamLUTs data frame.
#' @param beams the beams data frame.
#' @param x,y,z the coordinates of the 3D array.
#' @param variable name vector to override the default name of the variables (c('Alpha[Gy^(-1)]', 'Beta[Gy^(-2)]')).
#' @return A values object (if  x, y and z are specified), otherwise a sparse array.
#' @family BeamLUT
#' @export
get.alpha.beta.beamLUT <- function(beamLUTs, beams, x=NULL, y=NULL, z=NULL, variables=c('Alpha[Gy^(-1)]', 'Beta[Gy^(-2)]'))
{
  dose.bl <- beamLUTs$`DosePerEvent[Gy/primary]` * beams$fluence[beamLUTs$beamID]
  alpha.bl <- beamLUTs$`LethAlphaPerEvent[1/primary]` * beams$fluence[beamLUTs$beamID]
  beta.bl <- beamLUTs$`SqrtLethBetaPerEvent[1/primary]` * beams$fluence[beamLUTs$beamID]

  alpha.beta.bl <- aggregate(list(dose=dose.bl, alpha=alpha.bl, beta=beta.bl), by = list(voxelID=beamLUTs$voxelID), sum)


  if(is.null(x) | is.null(y) | is.null(z)) {
    alpha.beta.bl <- data.frame(alpha.beta.bl$voxelID, alpha.beta.bl$alpha/alpha.beta.bl$dose, alpha.beta.bl$beta^2/alpha.beta.bl$dose^2)
    names(alpha.beta.bl) <- c('voxelID', variables)
    return(alpha.beta.bl)
  } else {
    alpha <- beta <- array(NA, dim=c(length(x), length(y), length(z)))
    alpha[alpha.beta.bl$voxelID] = alpha.beta.bl$alpha/alpha.beta.bl$dose
    beta[alpha.beta.bl$voxelID] = alpha.beta.bl$beta^2/alpha.beta.bl$dose^2
    return(add.array.values(values=create.values(array.values = alpha, variables = variables[1], x = x, y = y, z = z), new.array=beta, variable=variables[2]))
  }
}

#' Get equivalent dose from beamLUTs (pure-dek)
#'
#' Get the net equivalent (biological) dose distribution (values object) from the beamLUTs.
#' @param beamLUTs the beamLUTs data frame.
#' @param beams the beams data frame.
#' @param x,y,z the coordinates of the 3D array.
#' @param alphaX,betaX the Linear-Quadratic parameters of the reference radiation (they can be two arrays/vectors)
#' @param variable to override the default name of the variable ('Dose[Gy]').
#' @return a values object.
#' @family BeamLUT
#' @export
get.bdose.beamLUT <- function(beamLUTs, beams, x, y, z, variable='BiologicalDose[Gy(RBE)]', alphaX, betaX)
{
  rbe <- get.rbe.beamLUT(beamLUTs = beamLUTs, beams = beams, x = x, y = y, z = z, alphaX = alphaX, betaX = betaX)
  dose <- get.dose.beamLUT(beamLUTs = beamLUTs, beams = beams, x = x, y = y, z = z)

  return(create.values(array.values = rbe$values*dose$values, variables = variable, x = x, y = y, z = z))
}

