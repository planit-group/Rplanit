# VOIS -------------------------------------------------------------------------


# Get the VOI index
#
#' Get the VOI index from the contours file
#' 
#' @param voi the VOI name.
#' @param contours a contours object (data frame).
#' @param file.contours the contours file.
#' @return the VOI index
#' 
#' @family VOIs
#' @export
get.voiindex <- function(voi, file.contours=NULL, contours=NULL)
{
  
  if(!is.null(contours)) {
    v <- unique(contours$id[contours$contour==voi])
    if(length(v)==0) {stop('selected voi not present in contours.')} else {return(v)}
  }
  else if(!is.null(file.contours)) {
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
  
}


#' Get VOIs from plan
#' 
#' Get the VOIs (\code{vois} object) from a \code{plan} object.
#' 
#' The \code{vois} object is similar (but different) to the values object.
#'
#' @param plan The \code{plan} object or a list of plan objects.
#' @param input Return the input vois (if exists).
#' @return the \code{vois} object (or a list of vois objects); The vois object is a list consisting of
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
get.vois <- function(plan, input=FALSE)
{
  # lista di piani
  if(class(plan)=='list') {
    vois <- list()
    for(i in 1:length(plan)) {
      vois[[i]] <- get.vois(plan[[i]])
    }
    return(vois)
  }
  
  # singolo piano
  else if(!is.null(plan[['vois']])) {
    return(plan[['vois']])
  } else {
    if(input) {vois.file <- plan[['inputVoisFile']]}
    else {vois.file <- plan[['outputVoisFile']]}
    if(is.null(vois.file)) {
      stop('vois object not present')
    } else {
      return(read.vois(vois.file))
    }
  }
}


#' Evaluate VOI logical matrix
#' 
#' Evaluate a boolean 3D array representing the volumetric occupancy of a specific VOI.
#' 
#' @param vois the \code{vois} object
#' @param voi the name of the VOI
#' @param Nv number of values (optional). This is used to produce a 4D array to account to variable multiplicity in the voxels (to be used directly with values objects with multiple variables)
#' @return boolean array representing the volumetric occupancy of a specific VOI.
#' 
#' @family VOIs
#' @export
get.voi.logical <- function(vois, voi, Nv=NULL) {
  # indice VOI
  #if(is.null(voi)) {cat('No VOI...\n'); return(0)}
  v <- which(vois$vois==voi)
  if(is.null(Nv)) {
    return(array(substr(vois$values, v, v)=='1', dim=c(vois$Nx, vois$Ny, vois$Nz)))
  } else {
    return(array(rep(substr(vois$values, v, v)=='1', each=Nv), dim=c(Nv, vois$Nx, vois$Ny, vois$Nz)))
  }
}

#' Get a subset of a  \code{vois} object.
#' @param values the \code{vois} object
#' @param xlim,ylim,zlim The ranges of the subset coordinates.
#' @family VOIs
#' @export
get.subset.vois <- function(vois, xlim=c(-Inf,Inf), ylim=c(-Inf,Inf), zlim=c(-Inf,Inf))
{
  # estrai range
  xlim[1] <- max(min(vois$x), xlim[1]); xlim[2] <- min(max(vois$x), xlim[2])
  ylim[1] <- max(min(vois$y), ylim[1]); ylim[2] <- min(max(vois$y), ylim[2])
  zlim[1] <- max(min(vois$z), zlim[1]); zlim[2] <- min(max(vois$z), zlim[2])
  
  ix.min <- max(which(vois$x<xlim[1]), 1); ix.max <- min(which(vois$x>xlim[2]), vois$Nx); # print(c(ix.min,ix.max))
  iy.min <- max(which(vois$y<ylim[1]), 1); iy.max <- min(which(vois$y>ylim[2]), vois$Ny)
  iz.min <- max(which(vois$z<zlim[1]), 1); iz.max <- min(which(vois$z>zlim[2]), vois$Nz)
  
  x <- vois$x[ix.min:ix.max]
  y <- vois$x[iy.min:iy.max]
  z <- vois$z[iz.min:iz.max]
  

  array.v <- vois$values[ix.min:ix.max,iy.min:iy.max,iz.min:iz.max]
  
  return(list(values = array.v, vois = vois$vois, x=x, y=y, z=z, Nx=length(x), Ny=length(y), Nz=length(z)))
}

#' Check if point(s) is(are) inside a VOI
#' @param x,y,z coordinate(s) of the point(s).
#' @param contours Contours object.
#' @param voi The VOI name.
#' @return A logical vector.
#' @family VOIs
#' @export
is.in.voi <- function(x, y, z, contours, voi, progressbar=FALSE)
{
  # seleziona voi ed identifica z
  zs.all <- sort(unique(contours$z)); dz <- mean(diff(zs.all)); zlim <- range(zs.all) + c(-dz,dz)/2
  contours <- subset(contours, contour==voi)
  if(nrow(contours)==0) {stop('voi ', voi, ' not present in contours')}
  zs <- sort(unique(contours$z))
  xlim <- range(contours$x)
  ylim <- range(contours$y)
  
  # rimuove punti out-of-box
  index.in <- (x>xlim[1]) & (x<xlim[2]) & (y>ylim[1]) & (y<ylim[2]) & (z>zlim[1]) & (z<zlim[2])
  x.in <- x[index.in]
  y.in <- y[index.in]
  z.in <- z[index.in]
  

  Np <- length(x.in)
  if(Np==0) return(rep(FALSE, length(x)))
  #print(xlim); print(ylim); print(zlim)
  message('evaluating ', Np, ' points of ', length(x), '...')
  inside <- rep(FALSE, Np)
  if(progressbar) {pb <- txtProgressBar(min = 0, max = Np, style = 3)}
  for(ip in 1:Np) {
    if(progressbar) {setTxtProgressBar(pb, ip)}
    
    # check out-of-box
    #if((x[ip]<xlim[1]) | (x[ip]>xlim[2]) | (y[ip]<ylim[1]) | (y[ip]>ylim[2]) | (z[ip]>zlim[2]) | (z[ip]<zlim[1])) {next}
    
    # check su z
    zp <- zs.all[ abs(z.in[ip]-zs.all) ==  min(abs(z.in[ip]-zs.all)) ]
    if( sum((zp %in% zs))==0 ) {next}
    
    # check su x,y
    zp <- zp[zp %in% zs][1]
    contours.z <- subset(contours, z==zp)
    x1 <- contours.z$x; x2 <- c(contours.z$x[2:(length(contours.z$x))], contours.z$x[1]) # segmenti
    y1 <- contours.z$y; y2 <- c(contours.z$y[2:(length(contours.z$y))], contours.z$y[1])
    #print(x1);print(x2)
    #print(y1);print(y2)
    
    # elimina segmenti verticali
    index.x.ok <- x1 != x2
    x1 <- x1[index.x.ok]
    x2 <- x2[index.x.ok]
    y1 <- y1[index.x.ok]
    y2 <- y2[index.x.ok]
    
    # identifica y0 sui segmenti
    y0 <- (y2-y1)*(x1-x.in[ip])/(x1-x2) + y1
    intersezioni <- (y.in[ip]<=y0) &
      ( (x.in[ip]>=x1 & x.in[ip]<=x2) | (x.in[ip]<=x1 & x.in[ip]>=x2))
    
    if( (sum(intersezioni) %% 2)>0 ) {inside[ip] <- TRUE}
    #if( sum(intersezioni)>1 ) {    print(intersezioni);stop()}
  }
  if(progressbar) {close(pb)}
  inside.all <- rep(FALSE, length(x))
  inside.all[index.in] <- inside
  return(inside.all)
}


#' Create a VOIS object from a contours object
#' @param contours The contours object.
#' @param x,y,z The coordinates of the x,y,z axes for the VOIS 3D grid.
#' @param vois Vois names vector. Optionally it is possible to specificy a specific set of vois to use.
#' @return a VOIS object.
#' @family VOIs
#' @export
create.vois <- function(contours, x, y, z, vois=NULL) {

  if(is.null(vois)) {vois <- unique(contours$contour)}
  else {vois <- vois[which(vois %in% unique(contours$contour))]}
  xyz <- expand.grid(list(x=x, y=y, z=z))
  
  Nv <- length(vois)
  for(iv in 1:Nv) {
    message('creating vois for ', vois[iv], ' ...')
    voi.logical <- is.in.voi(x = xyz$x, y = xyz$y, z = xyz$z, contours = contours, voi = vois[iv], progressbar = TRUE)
    voi.logical.char <- rep('0', nrow(xyz))
    voi.logical.char[voi.logical] <- '1'
    if(iv==1) {
      vois.logical.char.all <- voi.logical.char
    } else {
      vois.logical.char.all <- paste0(vois.logical.char.all, voi.logical.char)
    }
  }
  # dimensioni giuste
  dim(vois.logical.char.all) <- c(length(x), length(y), length(z))
  
  return(list(values=vois.logical.char.all, vois=vois, x=x, y=y, z=z, Nx=length(x), Ny=length(y), Nz=length(z)))
}

#' Create a VOIS object from a logical array
#' @param logical.array The logical array (TRUE -> the voi is present in the voxel.
#' @param x,y,z The coordinates of the x,y,z axes for the VOIS 3D grid.
#' @param voi the VOI name.
#' @return a VOIS object.
#' @family VOIs
#' @export
create.vois.from.logical.array <- function(logical.array, x, y, z, voi) {
  
  xyz <- expand.grid(list(x=x, y=y, z=z))
  
  message('creating vois for ', voi, ' ...')
  vois.logical.char <- rep('0', nrow(xyz))
  vois.logical.char[logical.array] <- '1'
  
  # dimensioni giuste
  dim(vois.logical.char) <- c(length(x), length(y), length(z))
  
  return(list(values=vois.logical.char, vois=voi, x=x, y=y, z=z, Nx=length(x), Ny=length(y), Nz=length(z)))
  
}

#' Combine 2 existing VOIS object into one
#' @param vois,vois.new the VOIS objects to be combined.
#' @return a VOIS object.
#' @family VOIs
#' @export
add.vois <- function(vois, vois.new)
{
  if(!any(vois$x==vois.new$x) | !any(vois$y==vois.new$y) | !any(vois$z==vois.new$z)) {
    stop('vois dimensions not congruent.')
  }
  if(any(vois$vois == vois.new$vois)) {
    stop('Some vois already present.')
  }
  
  vois.logical.char <- vois$values
  vois.logical.char.new <- vois.new$values
  vois.logical.char.all <- paste0(vois.logical.char, vois.logical.char.new)
  dim(vois.logical.char) <- c(length(vois$x), length(vois$y), length(vois$z))
  
  return(list(values=vois.logical.char.all, vois=c(vois$vois, vois.new$vois), x=vois$x, y=vois$y, z=vois$z, Nx=vois$Nx, Ny=vois$Ny, Nz=vois$Nz))
}

# CONTOURS ---------------------------------------------------------------------

#' Get contours from plan
#' 
#' Get the contours dataframe from the \code{plan} object.
#' 
#' @param plan the \code{plan} object. It can be a list of plans.
#' @return A \code{contours} dataframe (or a list of dataframes) consisting of:
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

#' @family Contours
#' @export
get.contours.list <- function(plans) {
  contours <- list()
  for(i in 1:length(plans)) {
    contours[[i]] <- get.contours(plans[[i]])
  }
  return(contours)
}

#' @family Contours
#' @export
get.contours.plankit.plan <- function(plan) {
  
  # prende direttamente data.frame contours se presente in plan
  if(!is.null(plan[['contours']])) return(plan[['contours']])
  
  # recupera informazioni da plan
  file.CT <- plan$ctFile
  file.contours <- plan$contoursFile
  targets <- plan$fields$targetVOIIndex
  
  return(read.contours(file.contours = file.contours, file.CT = file.CT))
  
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
#' @param voi the voi(s). If it is not specified, the volume of all VOIs is evaluated (optional)
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

#' Correct consistency of contours
#' 
#' Check and correct the consistency of a contours object (such as not uniform and/or ordered polygon numbering, negative index, etc.)
#' 
#' @param contours the contours dataframe
#' 
#' @export
#' @family Contours
sanitize.contours <- function(contours) {
  
  # uniformizza id
  ids <- sort(unique(contours$id))
  contours$id <- sapply(contours$id, function(v) which(v==ids)-1)
  
  # ordina
  contours <- contours[with(contours, order(id)), ]
  
  # ciclo sui voi
  for(iid in 1:length(ids)) {
    contour.name <- unique(contours$contour[contours$id==(iid-1)])
    message('sanitizing contour: ', contour.name)
    contours.voi <- subset(contours, contour==contour.name)
    
    # uniformizza poligoni
    polygons <- sort(unique(contours.voi$polygon))
    contours.voi$polygon <- sapply(contours.voi$polygon, function(v) which(v==polygons)-1)
    
    if(iid==1) {cc <- contours.voi} else {cc <- rbind(cc, contours.voi)}
  }
  
  # check corrispondenza z <-> slice rispetto alla CT (DA FARE!)
  
  return(cc)
}


#' Add color definitions to contours
#' 
#' Colours will be added using the "raimbow" colour selection and will be stored in a new "display.colour" column of the dataframe.
#' 
#' @param contours the contours dataframe.
#' @return contours dataframe with colour definition for each contour.
#' 
#' @export
#' @family Contours
add.colours.contours <- function(contours)
{
  contorni <- as.character(unique(contours$contour))
  colori <- rainbow(length(contorni))
  names(colori) <- contorni
  contours$display.color <- colori[as.character(contours$contour)]

  return(contours)
}


# R/W CONTOURS -----------------------------------------------------------------


#' Read contours from file
#' 
#' Read the contours from a "*.contours" file format (PlanKIT format) and return a \code{contours} dataframe.
#' 
#' @param file.contours the name of the contours file
#' @param file.CT the name of the CT file. The CT file (in PlanKIT format) is needed to properly reconstruc the 3D coordinates of the contours (alternatively to CT and z.CT)
#' @param CT the CT object (alternatively to file.CT and z.CT)
#' @param z.CT a vector for the z coordinates of the slices of the reference CT (alternatively to file.CT and CT)
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
read.contours <- function(file.contours, file.CT=NULL, CT=NULL, z.CT=NULL) {
  
  z <- NULL
  
  if(!is.null(z.CT)) {
    message('using z from z.CT...')
    z <- z.CT
  } else if(!is.null(CT)) {
    message('using z from CT...')
    z <- CT$z
  } else if(!is.null(file.CT)){
    message('using z from ', file.CT, '...')
    
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
  }
  
  if(is.null(z)) {stop('No z coordinates for contours')}
  
  
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
  
  # check per inconsistenze nella tabella
  i.wrong <- which(contours.ct$slice <0) # indici di punti "strani"
  Nwrong <- length(i.wrong)
  if(Nwrong>0) {
    warning('warning: unusual data in contours file... removing them.')
    print(summary(contours.ct))
    contours.ct <- contours.ct[-i.wrong,]
  }
  
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
  contours.ct$type <- contours.ct$type
  
  # trasforma fattori in strings
  contours.ct$contour <- as.character(contours.ct$contour)
  contours.ct$tissue <- as.character(contours.ct$tissue)
  
  
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
              file=paste0(name, '.contours'),
              col.names=FALSE, row.names=FALSE, append=TRUE)
  
  message('contours saved in file ', paste0(name, '.contours'))
}


#' Read contours from DICOM file
#' 
#' Read the contours from a DICOM file and return a \code{contours} data frame. The CT data is also mandatory
#' to identify the slice for each contour polygon. The CT data can be specified with one of the following:
#' folder.CT.dicom, CT or z.CT. z.CT has precendence to CT that has precendence to folder.CT.dicom.
#' 
#' @param file.contours.dicom the name of the contours file
#' @param folder.CT.dicom the name of the DICOM CT folder. The CT is needed to properly identify the slice for each contour polygon.
#' @param CT the CT object (alternative to folder.CT.dicom and z.CT)
#' @param z.CT a vector for the z coordinates of the slices of the reference CT (alternative to folder.CT.dicom and CT)
#' @return A \code{contours} dataframe consisting of:
#' \item{id}{index of the contours}
#' \item{polygon}{index of the polygon}
#' \item{slice}{index of the slice}
#' \item{x,y,z}{x,y,z coordinates}
#' \item{contour}{contour (VOI) name}
#' \item{tissue}{tissue+biological model name}
#' \item{type}{contour type (target, OAR, etc.)}
#' \item{colour}{the hex specification of the colour associated to the contour.}
#' 
#' @family R/W Contours
#' @export
read.contours.dicom <- function(file.contours.dicom, folder.CT.dicom=NULL, CT=NULL, z.CT=NULL, tissue='R3PIDEd0.2v0+MKM') {
  
  # CT data
  x_verspr <- y_versor <- 1
  if(is.null(z.CT)) {
    if(is.null(CT)) {
      CT <- read.3d.dicom(folder.CT.dicom, recursive=FALSE)
    }
    z.CT <- CT$z
  }
  
  
  # SS data
  header <- header2list(file = file.contours.dicom)
  
  message('creating contours data frame...')
  
  structureSetItems <- names(header[['StructureSetROISequence']])
  RTROIObservationsItems <- names(header[['RTROIObservationsSequence']])
  if (length(structureSetItems) != length(RTROIObservationsItems)) stop('VOI and observations not equal!')
  
  VOINumber <- VOIName <- VOIInterpretation <- contours <- NULL
  
  for(i in 1:length(structureSetItems)) {
    structureSetItem <- header[['StructureSetROISequence']][[structureSetItems[i]]]
    VOINumber[i] <- structureSetItem[['ROINumber']]
    VOIName[i] <- structureSetItem[['ROIName']]
    VOIInterpretation[i] = header[['RTROIObservationsSequence']][[RTROIObservationsItems[i]]][['RTROIInterpretedType']]
  }
  
  ROIContourItems <- names(header[['ROIContourSequence']])
  validContours <- rep(TRUE, length(ROIContourItems))
  
  for(n in 1:length(ROIContourItems)) {
    ROIContourItem <- header[['ROIContourSequence']][[ROIContourItems[n]]]
    refVOINumber = ROIContourItem[['ReferencedROINumber']]
    voiIndex <- which(VOINumber == refVOINumber)
    if( length(voiIndex) == 0 ) stop('contour index not present among VOI indices')
    
    contourItems <- names(ROIContourItem[['ContourSequence']]);
    for (p in 1:length(contourItems)) {
      contourItem <- ROIContourItem[['ContourSequence']][[contourItems[p]]]
      if ((contourItem[['ContourGeometricType']] != 'CLOSED_PLANAR') & (contourItem[['ContourGeometricType']] != 'POINT')) stop('Non closed planar contour')
      
      coordinates <- t(array(contourItem[['ContourData']], dim = c(3,length(contourItem[['ContourData']])/3)))
      z <- unique(coordinates[,3]) # deve essere uno solo in quanto "CLOSED_PLANAR"
      if(length(z)>1) stop(VOIName[voiIndex], ' not planar contour: z = ', z) # ma se il piano non è ortogonale a z...
      slice <- which(!(z <= z.CT | z > z.CT))
      if (length(slice) == 0) slice <- NA
      RGB <- ROIContourItem[['ROIDisplayColor']]
      colour <- rgb(RGB[1], RGB[2], RGB[3], maxColorValue = 255)
      
      contours.tmp <- data.frame(id=voiIndex, polygon=p, slice=slice,
                                 x=coordinates[,1], y=coordinates[,2], z=coordinates[,3],
                                 contour=VOIName[voiIndex], tissue=tissue, type=VOIInterpretation[voiIndex],
                                 colour=colour)
      contours <- rbind(contours, contours.tmp)
    }
    
  }
  
  return(contours)
}

read.3d.dicom <- function(dicom.folder, exclude=NULL, recursive=TRUE, verbose=TRUE, invert=TRUE, variable='HounsfieldNumber')
{
  dcmImages <- readDICOM(path=dicom.folder, verbose=verbose, recursive=recursive, exclude=exclude)
  dcm.info <- dicomTable(dcmImages$hdr)
  
  # trasforma le immagini...
  message('scaling images...')
  rs <- as.numeric(dcm.info$`0028-1053-RescaleSlope`)
  ri <- as.numeric(dcm.info$`0028-1052-RescaleIntercept`)
  for(i in 1:length(rs)) {
    dcmImages$img[[i]] <- dcmImages$img[[i]]*rs[i] + ri[i]
  }
  
  # inverte immagini
  if(invert) {
    message('inverting images...')
    for(i in 1:length(dcmImages$img)) {
      dcmImages$img[[i]] <- dcmImages$img[[i]][seq(dim(dcmImages$img[[i]])[1], 1), ]
    }
  }
  
  # ordina le immagini
  message('sorting images...')
  zz <- as.numeric(dcm.info$`0020-1041-SliceLocation`)
  oz <- order(zz)
  imgs <- list()
  for(i in 1:length(oz)) {
    imgs[[i]] <- dcmImages$img[[oz[i]]]
  }
  dcmImages$img <- imgs
  
  # array 3d
  Values.3d <- create3D(dcmImages); #return(Values.3d)
  Nx <- dim(Values.3d)[1]
  Ny <- dim(Values.3d)[2]
  Nz <- dim(Values.3d)[3]
  Nv <- 1
  
  # coordinates
  
  dxy <- as.numeric(unlist(strsplit(paste(dcm.info$`0028-0030-PixelSpacing`, collapse=' '), split=' ', fixed=TRUE)))
  
  vo <- as.numeric(unlist(strsplit(paste(dcm.info$`0020-0032-ImagePositionPatient`, collapse=' '), split=' ', fixed=TRUE)))
  voxel.origin <- c(mean(vo[seq(1, length(dxy)-2, by=3)]),
                    mean(vo[seq(2, length(dxy)-1, by=3)]),
                    mean(vo[seq(3, length(dxy), by=3)]))
  message('origin: ', voxel.origin[1], ', ', voxel.origin[2], ', ', voxel.origin[3])
  dx <- mean(dxy[seq(1, length(dxy)-1, by=2)])
  dy <- mean(dxy[seq(2, length(dxy), by=2)])
  x <- seq(from=voxel.origin[1], by=dx, length.out=Nx)
  y <- seq(from=voxel.origin[2], by=dy, length.out=Ny)
  z <- sort(zz)-min(zz) + voxel.origin[3]
  z <- sort(zz)
  
  # crea oggetto values
  values <- list(values=Values.3d, x=x, y=y, z=z, Nx=Nx, Ny=Ny, Nz=Nz, Nv=Nv, variables=variable)
  class(values) <- 'values'
  return(values)
}


# R/W VOIS ---------------------------------------------------------------------


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
read.vois <- function(file.name,
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

