# LUT --------------------------------------------------------------------------


#' Legge LUT
#'
#' Legge pencil-LUT da un file .4d e restituisce un oggetto lut (opzionalmente un dataframe)
#' @family LUT
#' @export
read.lut <- function(lut.name, dataframe=FALSE)
{

  library(reshape)
  
  lut.file <- paste(lut.name, '/', lut.name, '.lut.4d', sep='')
  #lut.file <- lut.name
  
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


#' Legge file zBP vs. E e ritorna una funzione d'interpolazione
#' 
#' @family LUT
#' @export
get.zBP.fun <- function(lut.name, EBP=FALSE)
{
  
  zBP.file <- paste(lut.name, '/', lut.name, '.zBraggPeaks.1d', sep='')
  
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


#' aggiungi alpha alla lut
#' 
#' la formula: LethAlphaPerEvent / DosePerEvent
#' 
#' @family LUT
#' @export
lut.add.alpha <- function(lut)
{
  message('evaluating alpha... ')
  if(class(lut)=='data.frame'){
    a <- subset(lut, variable=='LethAlphaPerEvent[1/primary]')
    d <- subset(lut, variable=='DosePerEvent[Gy/primary]')
    a$value <- a$value / d$value
    a$variable <- 'Alpha[Gy^(-1)]'
    return(rbind(lut, a))
  } else {
    Na <- which(lut$variables=='LethAlphaPerEvent[1/primary]')
    Nd <- which(lut$variables=='DosePerEvent[Gy/primary]')
    a <- lut$values[Na,,,,, drop=FALSE]
    d <- lut$values[Nd,,,,, drop=FALSE]
    a <- a / d
    #plot(a)
    av <- aperm(lut$values, c(2,3,4,5,1))
    a <- aperm(a, c(2,3,4,5,1))
    av <- array(c(av,a), c(lut$NE, lut$Nx, lut$Ny, lut$Nz, lut$Nv+1) )
    lut$values <- aperm(av, c(5,1,2,3,4))
    lut$variables[lut$Nv+1] <- 'Alpha[Gy^(-1)]'
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

#' Evaluate BeamLUT 
#' 
#' Evaluate the beamLUT of the selected beam(s) for the specific plan.
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


#' scrive le matrici V_{b,v} (beamLUT) 
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
    df.s <- get.sparse.array.from.values(values=vb, variable=variable, threshold=threshold)
    #df.s$value <- df.s$value/beams$fluence[b]
    df.s$beamID <- beams$beamID[b]
    
      if(b==1) {append=FALSE} else {append=TRUE}
      write.table(df.s, file=file.name, col.names=FALSE, row.names=FALSE, append=append)
    
  }
}

