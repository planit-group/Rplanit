# VOIS -------------------------------------------------------------------------


# Get the VOI index
#
#' Get the VOI index from the contours file
#' 
#' @param voi the VOI name
#' @param file.contours the contours file
#' @return the VOI index
#' 
#' @family VOIs
#' @export
get.voiindex <- function(voi, file.contours=NULL)
{
  
  # header
  file.con <- file(file.contours, "rb") # read binary
  myline <- readLines(file.con, n=1)
  myline.splitted <- unlist(strsplit(myline[1], ' +'))
  Nc <- as.numeric(myline.splitted[1])
  con.names <- myline.splitted[2:length(myline.splitted)]
  close(file.con)
  
  # splitta tessuto+modello
  con <- tis <- con.names
  Nc <- length(con.names)
  con.name <- strsplit(con.names, '_', fixed=TRUE)
  for (i in 1:Nc) {
    con[i] <- unlist(con.name[i])[1]
    tis[i] <- unlist(con.name[i])[2]
  }
  
  
  v <- which(con==voi)
  return(v-1)
  
}


#' Get VOIs from plan
#' 
#' Get the VOIs (\code{vois} object) from a \code{plan} object.
#' 
#' The \code{vois} object is similar (but different) to the values object.
#'
#' @param plan the \code{plan} object
#' @return the \code{vois} object; a list consisting of
#' \item{values}{3D "character" array (logical strings such as "10010110")}
#' \item{vois}{vector of VOI names}
#' \item{x}{vector of x coordinates}
#' \item{y}{vector of y coordinates}
#' \item{z}{vector of z coordinates}
#' \item{Nx}{number of voxels along x}
#' \item{Ny}{number of voxels along y}
#' \item{Nz}{number of voxels along z}
#' \item{Nv}{number of VOIs}
#' \item{file}{reference file}
#' 
#' @family VOIs
#' @export
get.vois <- function(plan)
{
  if(is.null(plan[['voisFile']])) {
    cat('The plan "', plan[['name']], '" has no vois file.\n', sep='')
    return(NULL)
  } 
  return(read.vois.array( paste(plan[['name']], '/', plan[['voisFile']], sep='') ))
}


#' Evaluate VOI logical matrix
#' 
#' Evaluate a boolean 3D array representing the volumetric occupancy of a specific VOI.
#' 
#' @param vois the \code{vois} object
#' @param voi the name of the VOI
#' @return boolean array representing the volumetric occupancy of a specific VOI.
#' 
#' @family VOIs
#' @export
get.voi.logical <- function(vois, voi=NULL) {
  # indice VOI
  if(is.null(voi)) {cat('No VOI...\n'); return(0)}
  v <- which(vois$vois==voi)
  return <- substr(vois$values, v, v)=='1'
}


# CONTOURS ---------------------------------------------------------------------

#' Get contours from plan
#' 
#' Get the contours dataframe from the \code{plan} object.
#' 
#' @param plan the \code{plan} object
#' @return A \code{contours} dataframe consisting of:
#' \item{id}{index of the contours}
#' \item{polygon}{index of the polygon}
#' \item{slice}{index of the slice}
#' \item{x,y,z}{x,y,z coordinates}
#' \item{contour}{contour (VOI) name}
#' \item{tissue}{tissue+biological model name}
#' \item{type}{contour type (target, OAR, etc.)}
#' 
#' @family Contours
#' @export
get.contours <- function(plan) UseMethod('get.contours')

get.contours.plankit.plan <- function(plan) {
  
  # recupera informazioni da plan
  file.CT <- plan$ctFile
  file.contours <- plan$contoursFile
  targets <- plan$fields$targetVOIIndex
  
  
  cat('reading contours: ', file.contours, ', ', file.CT, '\n', sep='')
  
  # apri connsessione
  file.3d <- file(file.CT, "rb") # read binary
  
  # leggi l'header CT
  cat('reading header CT...\n')
  myline <- readLines(file.3d, n=8) # legge le prime 8 linee, ogni linea e' un elemento del vettore
  close(file.3d)
  
  # parsing dell'header
  # splitta la stringa in sottostringhe delimitate da uno o piu' spazi (regexpr: " +")
  myline.splitted <- unlist(strsplit(myline[1], ' +'))
  Nx <- as.numeric(myline.splitted[1])
  myline.splitted <- unlist(strsplit(myline[2], ' +'))
  x <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[3], ' +'))
  Ny <- as.numeric(myline.splitted[1])
  myline.splitted <- unlist(strsplit(myline[4], ' +'))
  y <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[5], ' +'))
  Nz <- as.numeric(myline.splitted[1])
  myline.splitted <- unlist(strsplit(myline[6], ' +'))
  z <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[7], ' +'))
  Nv <- as.numeric(myline.splitted[1])
  values <- myline.splitted[2:length(myline.splitted)]
  
  Ntot <- Nx * Ny * Nz * Nv
  
  # trasforma intervalli in coordinate puntuali
  x <- (x[1:Nx] + x[2:(Nx+1)])/2
  y <- (y[1:Ny] + y[2:(Ny+1)])/2
  z <- (z[1:Nz] + z[2:(Nz+1)])/2
  
  
  # legge contorni
  
  # header
  file.con <- file(file.contours, "rb") # read binary
  myline <- readLines(file.con, n=1)
  myline.splitted <- unlist(strsplit(myline[1], ' +'))
  Nc <- as.numeric(myline.splitted[1])
  con.names <- myline.splitted[2:length(myline.splitted)]
  close(file.con)
  
  print(con.names)
  
  # body
  contours.ct <- read.table(file.contours, skip=1)
  names(contours.ct) <- c('id', 'polygon', 'slice', 'x', 'y')
  contours.ct$z <- (sapply(contours.ct$slice, function(v) z[v+1]))
  contours.ct$contour <- factor(sapply(contours.ct$id, function(v) con.names[v+1]), levels=con.names)
  
  # identifica tessuto+modello
  con.names <- levels(contours.ct$contour)
  con <- tis <- con.names
  Nc <- length(con.names)
  con.name <- strsplit(con.names, '_', fixed=TRUE)
  for (i in 1:Nc) {
    con[i] <- unlist(con.name[i])[1]
    tis[i] <- unlist(con.name[i])[2]
  }
  levels(contours.ct$contour) <- con
  contours.ct$tissue <- contours.ct$contour
  levels(contours.ct$tissue) <- tis
  
  # identifica target
  contours.ct$type <- 'VOI/OAR'
  contours.ct$type[contours.ct$id %in% targets] <- 'PTV'
  contours.ct$type <- as.factor(contours.ct$type)
  
  return(contours.ct)
  
}

get.contours.gate.plan <- function(plan) {
  if(is.null(plan$contours)) {
    message('Contours not present in Gate plan')
    return(NULL)
  }
  return(plan$contours)
}


#' Evaluate volume from contours
#' 
#' Evaluate the volume defined by its contours.
#' 
#' @param contours the contours dataframe
#' @voi the voi(s). If it is not specified, the volume of all VOIs is evaluated (optional)
#' 
#' @export
#' @importFrom splancs areapl
#' @family Contours
volume.contours <- function(contours, voi=NULL)
{
  
  # rimuovi eventuali NA
  contours <- subset(contours, !is.na(z))
  
  # selezione vois
  if(is.null(voi)) {voi <- unique(contours$contour)}
  
  N.voi <- length(voi)
  for(i in 1:N.voi) {
    cont.voi <- subset(contours, contour==voi[i])
    zs <- unique(cont.voi$z)
    dzs <- diff(points2intervals(zs))
    volume <- 0
    for(j in 1:length(zs)) {
      cont.voi.zs <- subset(cont.voi, z==zs[j])
      x <- cont.voi.zs$x
      y <- cont.voi.zs$y
      pol <- cbind(x,y)
      volume <- volume + areapl(pol)
    }
    volume.df.tmp <-  data.frame(contour=voi[i], volume=volume)
    if(i==1) {
      volume.df <- volume.df.tmp
    } else {
      volume.df <- rbind(volume.df, volume.df.tmp)
    }
  }
  
  return(volume.df) 
}

# R/W CONTOURS -----------------------------------------------------------------


#' Read contours from file
#' 
#' Read the contours from a "*.contours" file format (PlanKIT format) and return a \code{contours} dataframe.
#' 
#' @param file.contours the name of the contours file
#' @param file.CT the name of the CT file. The CT file (in PlanKIT format) is needed to properly reconstruc the 3D coordinates of the contours.
#' @return A \code{contours} dataframe consisting of:
#' \item{id}{index of the contours}
#' \item{polygon}{index of the polygon}
#' \item{slice}{index of the slice}
#' \item{x,y,z}{x,y,z coordinates}
#' \item{contour}{contour (VOI) name}
#' \item{tissue}{tissue+biological model name}
#' \item{type}{contour type (target, OAR, etc.)}
#' 
#' @family R/W Contours
#' @export
read.contours <- function(file.contours, file.CT) {
  
  # leggi header CT
  cat('reading contours: ', file.contours, ', ', file.CT, '\n', sep='')
  
  # apri connessione
  file.3d <- file(file.CT, "rb") # read binary
  
  # leggi l'header
  cat('reading header CT...\n')
  myline <- readLines(file.3d, n=8) # legge le prime 8 linee, ogni linea e' un elemento del vettore
  close(file.3d)
  
  # parsing dell'header
  # splitta la stringa in sottostringhe delimitate da uno o piu' spazi (regexpr: " +")
  myline.splitted <- unlist(strsplit(myline[1], ' +'))
  Nx <- as.numeric(myline.splitted[1])
  myline.splitted <- unlist(strsplit(myline[2], ' +'))
  x <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[3], ' +'))
  Ny <- as.numeric(myline.splitted[1])
  myline.splitted <- unlist(strsplit(myline[4], ' +'))
  y <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[5], ' +'))
  Nz <- as.numeric(myline.splitted[1])
  myline.splitted <- unlist(strsplit(myline[6], ' +'))
  z <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[7], ' +'))
  Nv <- as.numeric(myline.splitted[1])
  values <- myline.splitted[2:length(myline.splitted)]
  
  Ntot <- Nx * Ny * Nz * Nv
  cat('number of voxels:', Nv, 'x', Nx, 'x', Ny, 'x', Nz, '=', Ntot, '\n')
  cat('variables:', values, '\n')
  
  # trasforma intervalli in coordinate puntuali
  x <- (x[1:Nx] + x[2:(Nx+1)])/2
  y <- (y[1:Ny] + y[2:(Ny+1)])/2
  z <- (z[1:Nz] + z[2:(Nz+1)])/2
  
  
  # legge contorni
  
  cat('reading contours...\n')
  
  # header
  file.con <- file(file.contours, "rb") # read binary
  myline <- readLines(file.con, n=1)
  myline.splitted <- unlist(strsplit(myline[1], ' +'))
  Nc <- as.numeric(myline.splitted[1])
  con.names <- myline.splitted[2:length(myline.splitted)]
  close(file.con)
  
  print(con.names)
  
  # body
  contours.ct <- read.table(file.contours, skip=1)
  names(contours.ct) <- c('id', 'polygon', 'slice', 'x', 'y')
  contours.ct$z <- (sapply(contours.ct$slice, function(v) z[v+1]))
  contours.ct$contour <- factor(sapply(contours.ct$id, function(v) con.names[v+1]), levels=con.names)
  
  # identifica tessuto+modello
  con.names <- levels(contours.ct$contour)
  con <- tis <- con.names
  Nc <- length(con.names)
  con.name <- strsplit(con.names, '_', fixed=TRUE)
  for (i in 1:Nc) {
    con[i] <- unlist(con.name[i])[1]
    tis[i] <- unlist(con.name[i])[2]
  }
  levels(contours.ct$contour) <- con
  contours.ct$tissue <- contours.ct$contour
  levels(contours.ct$tissue) <- tis
  
  # identifica target
  contours.ct$type <- 'VOI/OAR'
  contours.ct$type[substr(contours.ct$contour,1,3) %in% 'PTV'] <- 'PTV'
  contours.ct$type <- as.factor(contours.ct$type)

  
  return(contours.ct)
  
}


#' Save contours
#' 
#' Save the \code{contours} dataframe in a file with PlanKIT format ("*.contours")
#' 
#' @param contours the \code{contours} dataframe
#' @param name the name of the file (the contours will be saved in \code{<name>.contours})
#' 
#' @family R/W Contours
#' @export
write.contours <- function(contours, name)
{
  numberOfvois <- length(unique(contours$id))
  vois <- unique(paste(contours$contour, contours$tissue, sep='_'))
  
  con <- file(paste(name, '.contours', sep=''), "w") # open for writing in text mode
  writeLines(paste(numberOfvois, " "), con=con, sep='')
  writeLines(paste(vois, collapse=' '), con=con, sep='')
  writeLines('\n', con=con, sep='')
  close(con)
  write.table(contours[c('id', 'polygon', 'slice', 'x', 'y')],
              file=paste(name, '.contours', sep=''),
              col.names=FALSE, row.names=FALSE, append=TRUE)
}


# R/W VOIS ---------------------------------------------------------------------


#' Legge i vois da un file .3d e restituisce un dataframe "esteso"
#' 
#' OBSOLETE
#' 
#' @family R/W VOIs
#' @export
read.vois <- function(file.name,
                      x.min=-Inf, x.max=+Inf,
                      y.min=-Inf, y.max=+Inf,
                      z.min=-Inf, z.max=+Inf) {
  
  #################
  # leggi file vois
  #################
  
  cat('Reading vois:', file.name, '\n')
  
  # apri connsessione
  file.3d <- file(file.name, "rb") # read binary
  
  # leggi l'header
  cat('Reading header...\n')
  myline <- readLines(file.3d, n=8) # legge le prime 8 linee, ogni linea e' un elemento del vettore
  
  # parsing dell'header
  # splitta la stringa in sottostringhe delimitate da uno o piu' spazi (regexpr: " +")
  myline.splitted <- unlist(strsplit(myline[1], ' +'))
  Nx <- as.numeric(myline.splitted[1])
  myline.splitted <- unlist(strsplit(myline[2], ' +'))
  x <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[3], ' +'))
  Ny <- as.numeric(myline.splitted[1])
  myline.splitted <- unlist(strsplit(myline[4], ' +'))
  y <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[5], ' +'))
  Nz <- as.numeric(myline.splitted[1])
  myline.splitted <- unlist(strsplit(myline[6], ' +'))
  z <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[7], ' +'))
  Nv <- as.numeric(myline.splitted[1])
  voi.names <- myline.splitted[2:length(myline.splitted)]
  
  Ntot <- Nx * Ny * Nz * Nv
  cat('Number of voxels:', Nv, 'x', Nx, 'x', Ny, 'x', Nz, '=', Ntot, '\n')
  cat('VOIs:', voi.names, '\n')
  
  # trasforma intervalli in coordinate puntuali
  x <- (x[1:Nx] + x[2:(Nx+1)])/2
  y <- (y[1:Ny] + y[2:(Ny+1)])/2
  z <- (z[1:Nz] + z[2:(Nz+1)])/2
  
  # crea e legge array
  cat('Reading VOIs data...\n')
  vois.3d <- array(readLines(file.3d), dim=c(Nx, Ny, Nz))
  
  # rimuove gli spazi
  vois.3d <- gsub(' ', '', vois.3d)
  
  # chiude file
  close(file.3d)
  
  
  #######################
  # selezione sottovolume
  #######################
  
  cat('Selecting subvolume...\n')
  
  # trova estremi nelle coordinate disponibili
  xx.min <- min(x[x>=x.min])
  yy.min <- min(y[y>=y.min])
  zz.min <- min(z[z>=z.min])
  xx.max <- max(x[x<=x.max])
  yy.max <- max(y[y<=y.max])
  zz.max <- max(z[z<=z.max])
  
  i.min <- which(x==xx.min)
  j.min <- which(y==yy.min)
  k.min <- which(z==zz.min)
  i.max <- which(x==xx.max)
  j.max <- which(y==yy.max)
  k.max <- which(z==zz.max)
  
  cat('subvolume: (', xx.min, ', ', xx.max, ') ',
      '(', yy.min, ', ', yy.max, ') ',
      '(', zz.min, ', ', zz.max, ')\n', sep='')
  cat('indexes: (', i.min, ', ', i.max, ') ',
      '(', j.min, ', ', j.max, ') ',
      '(', k.min, ', ', k.max, ')\n', sep='')
  
  # seleziona sottoarray
  cat('dimensioni:' ,dim(vois.3d), '\n')
  vois.3d <- vois.3d[i.min:i.max, j.min:j.max, k.min:k.max]
  
  # trim delle coordinate
  x <- x[x>=xx.min & x<=xx.max]
  y <- y[y>=yy.min & y<=yy.max]
  z <- z[z>=zz.min & z<=zz.max]
  
  
  ############################
  # crea dataframe strutturato
  ############################
  
  cat('Melting data into data.frame...\n')
  library(reshape)
  vois <- melt(vois.3d)
  rm(vois.3d)
  colnames(vois) <- c('x', 'y', 'z', 'value')
  vois$x <- (sapply(vois$x, function(v) x[v]))
  vois$y <- (sapply(vois$y, function(v) y[v]))
  vois$z <- (sapply(vois$z, function(v) z[v]))
  
  # trasforma in caratteri
  vois$value <- as.character(vois$value)
  
  # per testing...
  # return(vois)
  
  # estende il data.frame con denominazione esplicita dei VOIs
  for (v in 1:Nv) {
    cat('extending VOI:', voi.names[v], '...\n')    
    index <- substr(vois$value, v, v)=='1'
    vois.temp <- data.frame(x=vois$x[index], y=vois$y[index], z=vois$z[index], VOI=voi.names[v])
    if (v==1) {
      vois.ext <- vois.temp
    } else {
      vois.ext <- rbind(vois.ext, vois.temp)
    }
  }
  
  return(vois.ext)  
}


#' Legge i VOIS da un file .3d e restituisce un dataframe "esteso".
#' 
#' La creazione del dataframe e' fatta a "splices" per ottimizzare la memoria.
#' 
#' OBSOLETE
#' 
#' @family R/W VOIs
#' @export
read.vois.spliced <- function(file.name,
                              x.min=-Inf, x.max=+Inf,
                              y.min=-Inf, y.max=+Inf,
                              z.min=-Inf, z.max=+Inf,
                              Nsplices=8) {
  library(reshape)
  
  #################
  # leggi file vois
  #################
  
  cat('Reading vois:', file.name, '\n')
  
  # apri connsessione
  file.3d <- file(file.name, "rb") # read binary
  
  # leggi l'header
  cat('Reading header...\n')
  myline <- readLines(file.3d, n=8) # legge le prime 8 linee, ogni linea e' un elemento del vettore
  
  # parsing dell'header
  # splitta la stringa in sottostringhe delimitate da uno o piu' spazi (regexpr: " +")
  myline.splitted <- unlist(strsplit(myline[1], ' +'))
  Nx <- as.numeric(myline.splitted[1])
  myline.splitted <- unlist(strsplit(myline[2], ' +'))
  x <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[3], ' +'))
  Ny <- as.numeric(myline.splitted[1])
  myline.splitted <- unlist(strsplit(myline[4], ' +'))
  y <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[5], ' +'))
  Nz <- as.numeric(myline.splitted[1])
  myline.splitted <- unlist(strsplit(myline[6], ' +'))
  z <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[7], ' +'))
  Nv <- as.numeric(myline.splitted[1])
  voi.names <- myline.splitted[2:length(myline.splitted)]
  
  Ntot <- Nx * Ny * Nz * Nv
  cat('Number of voxels:', Nv, 'x', Nx, 'x', Ny, 'x', Nz, '=', Ntot, '\n')
  cat('VOIs:', voi.names, '\n')
  
  # trasforma intervalli in coordinate puntuali
  x <- (x[1:Nx] + x[2:(Nx+1)])/2
  y <- (y[1:Ny] + y[2:(Ny+1)])/2
  z <- (z[1:Nz] + z[2:(Nz+1)])/2
  
  # crea e legge array
  cat('Reading VOIs data...\n')
  vois.3d <- array(readLines(file.3d), dim=c(Nx, Ny, Nz))
  
  # rimuove gli spazi
  vois.3d <- gsub(' ', '', vois.3d)
  
  # chiude file
  close(file.3d)
  
  
  #######################
  # selezione sottovolume
  #######################
  
  cat('Selecting subvolume...\n')
  
  # trova estremi nelle coordinate disponibili
  xx.min <- min(x[x>=x.min])
  yy.min <- min(y[y>=y.min])
  zz.min <- min(z[z>=z.min])
  xx.max <- max(x[x<=x.max])
  yy.max <- max(y[y<=y.max])
  zz.max <- max(z[z<=z.max])
  
  i.min <- which(x==xx.min)
  j.min <- which(y==yy.min)
  k.min <- which(z==zz.min)
  i.max <- which(x==xx.max)
  j.max <- which(y==yy.max)
  k.max <- which(z==zz.max)
  
  cat('subvolume: (', xx.min, ', ', xx.max, ') ',
      '(', yy.min, ', ', yy.max, ') ',
      '(', zz.min, ', ', zz.max, ')\n', sep='')
  cat('indexes: (', i.min, ', ', i.max, ') ',
      '(', j.min, ', ', j.max, ') ',
      '(', k.min, ', ', k.max, ')\n', sep='')
  
  # seleziona sottoarray
  cat('dimensioni:' ,dim(vois.3d), '\n')
  vois.3d <- vois.3d[i.min:i.max, j.min:j.max, k.min:k.max]
  
  # trim delle coordinate
  x <- x[x>=xx.min & x<=xx.max]
  y <- y[y>=yy.min & y<=yy.max]
  z <- z[z>=zz.min & z<=zz.max]
  Nx <- length(x)
  Ny <- length(y)
  Nz <- length(z)
  
  
  ############################
  # crea dataframe strutturato
  ############################
  
  # splicing
  cat('splicing...\n')
  
  # calcola numero di slices per splice
  Ns <- floor(Nz/Nsplices)
  Ns <- rep(Ns, Nsplices)
  Ns[1] <- Ns[1] + (Nz-sum(Ns))
  Ns.cum <- c(0,cumsum(Ns))
  
  initial <- TRUE
  for(i in 1:Nsplices) {
    
    k.min <- Ns.cum[i]+1
    k.max <- Ns.cum[i+1]
    
    cat('splice ', i, ' (Ns = ', Ns[i], ' from k = ', k.min, ' to k = ', k.max, ')...\n', sep='')
    
    # seleziona splice
    vois.3d.s <-vois.3d[i.min:i.max, j.min:j.max, k.min:k.max]
    
    # melting
    cat('melting...\n')
    vois.s <- melt(vois.3d.s)
    #print(summary(vois.s))
    
    # managing
    vois.s$X3 <- vois.s$X3 + k.min-1
    colnames(vois.s) <- c('x', 'y', 'z', 'value')
    # trasforma in caratteri
    vois.s$value <- as.character(vois.s$value)
    
    # estende il data.frame con denominazione esplicita dei VOIs
    for (v in 1:Nv) {
      cat('extending VOI:', voi.names[v])    
      index <- substr(vois.s$value, v, v)=='1'
      Nvoxels <- length(vois.s$x[index])
      cat(' (found', Nvoxels, 'voxels) ...\n')
      if(Nvoxels>0) {
        vois.s.temp <- data.frame(x=vois.s$x[index], y=vois.s$y[index], z=vois.s$z[index], VOI=voi.names[v])
        vois.s.temp$x <- (sapply(vois.s.temp$x, function(v) x[v]))
        vois.s.temp$y <- (sapply(vois.s.temp$y, function(v) y[v]))
        vois.s.temp$z <- (sapply(vois.s.temp$z, function(v) z[v]))
        if (initial) {
          vois.ext <- vois.s.temp
          initial <- FALSE
        } else {
          vois.ext <- rbind(vois.ext, vois.s.temp)
        }
      }
    }
    
  }
  
  return(vois.ext)  
}


#' Read VOIs from file
#' 
#' Read the VOIS from a file (VOIs *.3d PlanKIT format) and return a \code{vois} object.
#' 
#' @param file.name the file name
#' @param x.min,x.max,y.min,y.max,z.min,z.max load only the volume inside the specified ranges (optinal)
#' @return The \code{vois} object
#' 
#' @family R/W VOIs
#' @export
read.vois.array <- function(file.name,
                            x.min=-Inf, x.max=+Inf,
                            y.min=-Inf, y.max=+Inf,
                            z.min=-Inf, z.max=+Inf)
{
  
  cat('Reading vois:', file.name, '\n')
  
  # apri connsessione
  file.3d <- file(file.name, "rb") # read binary
  
  # leggi l'header
  cat('Reading header...\n')
  myline <- readLines(file.3d, n=8) # legge le prime 8 linee, ogni linea e' un elemento del vettore
  
  # parsing dell'header
  # splitta la stringa in sottostringhe delimitate da uno o piu' spazi (regexpr: " +")
  myline.splitted <- unlist(strsplit(myline[1], ' +'))
  Nx <- as.numeric(myline.splitted[1])
  myline.splitted <- unlist(strsplit(myline[2], ' +'))
  x <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[3], ' +'))
  Ny <- as.numeric(myline.splitted[1])
  myline.splitted <- unlist(strsplit(myline[4], ' +'))
  y <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[5], ' +'))
  Nz <- as.numeric(myline.splitted[1])
  myline.splitted <- unlist(strsplit(myline[6], ' +'))
  z <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[7], ' +'))
  Nv <- as.numeric(myline.splitted[1])
  voi.names <- myline.splitted[2:length(myline.splitted)]
  
  Ntot <- Nx * Ny * Nz * Nv
  cat('Number of voxels:', Nv, 'x', Nx, 'x', Ny, 'x', Nz, '=', Ntot, '\n')
  cat('VOIs:', voi.names, '\n')
  
  # trasforma intervalli in coordinate puntuali
  x <- (x[1:Nx] + x[2:(Nx+1)])/2
  y <- (y[1:Ny] + y[2:(Ny+1)])/2
  z <- (z[1:Nz] + z[2:(Nz+1)])/2
  
  # crea e legge array
  cat('Reading VOIs data...\n')
  vois.3d <- array(readLines(file.3d), dim=c(Nx, Ny, Nz))
  
  # rimuove gli spazi
  vois.3d <- gsub(' ', '', vois.3d)
  
  # chiude file
  close(file.3d)
  
  
  # selezione sottovolume
  
  cat('Selecting subvolume...\n')
  
  # trova estremi nelle coordinate disponibili
  xx.min <- min(x[x>=x.min])
  yy.min <- min(y[y>=y.min])
  zz.min <- min(z[z>=z.min])
  xx.max <- max(x[x<=x.max])
  yy.max <- max(y[y<=y.max])
  zz.max <- max(z[z<=z.max])
  
  i.min <- which(x==xx.min)
  j.min <- which(y==yy.min)
  k.min <- which(z==zz.min)
  i.max <- which(x==xx.max)
  j.max <- which(y==yy.max)
  k.max <- which(z==zz.max)
  
  cat('subvolume: (', xx.min, ', ', xx.max, ') ',
      '(', yy.min, ', ', yy.max, ') ',
      '(', zz.min, ', ', zz.max, ')\n', sep='')
  cat('indexes: (', i.min, ', ', i.max, ') ',
      '(', j.min, ', ', j.max, ') ',
      '(', k.min, ', ', k.max, ')\n', sep='')
  
  # seleziona sottoarray
  cat('dimensioni:' ,dim(vois.3d), '\n')
  vois.3d <- vois.3d[i.min:i.max, j.min:j.max, k.min:k.max]
  
  # trim delle coordinate
  x <- x[x>=xx.min & x<=xx.max]
  y <- y[y>=yy.min & y<=yy.max]
  z <- z[z>=zz.min & z<=zz.max]
  Nx <- length(x)
  Ny <- length(y)
  Nz <- length(z)
  
  # crea struttura array
  return(list(values=vois.3d, vois=voi.names, x=x, y=y, z=z, Nx=Nx, Ny=Ny, Nz=Nz, Nv=Nv, file=file.name))
  
}

