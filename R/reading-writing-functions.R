# R/W ARRAYS -------------------------------------------------------------------

#' Read 3d array (PlanKIT)
#' 
#' Read file with format "3d" (PlanKIT) and return a \code{values} object.
#' 
#' @param file.name the file name
#' @param x.min,x.max,y.min,y.max,z.min,z.max load only the volume inside the specified ranges (optional)
#' 
#' @family R/W Arrays
#' @export
read.3d <- function(file.name,
                          x.min=-Inf, x.max=+Inf,
                          y.min=-Inf, y.max=+Inf,
                          z.min=-Inf, z.max=+Inf)
  {
  
  cat('reading values:', file.name, '\n')
  
  # apri connsessione
  file.3d <- file(file.name, "rb") # read binary
  
  # leggi l'header
  cat('reading header...\n')
  myline <- readLines(file.3d, n=8) # legge le prime 8 linee, ogni linea e' un elemento del vettore
  
  # parsing dell'header
  # splitta la stringa in sottostringhe delimitate da uno o piu' spazi (regexpr: " +")
  myline.splitted <- unlist(strsplit(myline[1], ' +'))
  Nx <- as.numeric(myline.splitted[1])
  x.type <- myline.splitted[3]
  myline.splitted <- unlist(strsplit(myline[2], ' +'))
  x <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[3], ' +'))
  Ny <- as.numeric(myline.splitted[1])
  y.type <- myline.splitted[3]
  myline.splitted <- unlist(strsplit(myline[4], ' +'))
  y <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[5], ' +'))
  Nz <- as.numeric(myline.splitted[1])
  z.type <- myline.splitted[3]
  myline.splitted <- unlist(strsplit(myline[6], ' +'))
  z <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[7], ' +'))
  Nv <- as.numeric(myline.splitted[1])
  variables <- myline.splitted[2:length(myline.splitted)]
  
  type <- myline[8]
  
  Ntot <- Nx * Ny * Nz * Nv
  cat('number of voxels:', Nv, 'x', Nx, 'x', Ny, 'x', Nz, '=', Ntot, '\n')
  cat('variables:', variables, '\n')
  cat('type:', type, '\n')
  
  # trasforma intervalli in coordinate puntuali
  if(x.type=='INTERVAL') {x <- (x[1:Nx] + x[2:(Nx+1)])/2}
  if(y.type=='INTERVAL') {y <- (y[1:Ny] + y[2:(Ny+1)])/2}
  if(z.type=='INTERVAL') {z <- (z[1:Nz] + z[2:(Nz+1)])/2}
  
  
  # crea e legge array (4d)
  if(type=='BINARY') {
    cat('reading binary data...\n')
    Values.3d <- array(readBin(file.3d, numeric(), Ntot), dim=c(Nv, Nx, Ny, Nz))
  } else {
    cat('Reading value data...\n')
    Values.3d <- (array(as.numeric(unlist(strsplit(readLines(file.3d), ' +'))), dim=c(Nv, Nx, Ny, Nz)))
    #print(values)
    #array.1d <- array(readLines(file.1d), dim=c(Nv, Nx))
  }
  
  # chiudi connessione
  close(file.3d)
  
  #######################
  # selezione sottovolume
  #######################
  
  cat('selecting subvolume...\n')
  
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
  #Values.3d <- Values.3d[ 1:Nv, i.min:i.max, j.min:j.max, k.min:k.max]
  Values.3d <- Values.3d[ 1:Nv, i.min:i.max, j.min:j.max, k.min:k.max, drop=FALSE]
  if(Nv==1) {Values.3d <- array(Values.3d, dim=c(Nx,Ny,Nz))}
  
  # trim delle coordinate
  x <- x[x>=xx.min & x<=xx.max]
  y <- y[y>=yy.min & y<=yy.max]
  z <- z[z>=zz.min & z<=zz.max]
  Nx <- length(x)
  Ny <- length(y)
  Nz <- length(z)
  
  # crea oggetto values
  values <- list(values=Values.3d, x=x, y=y, z=z, Nx=Nx, Ny=Ny, Nz=Nz, Nv=Nv, variables=variables, file=file.name)
  class(values) <- 'values'
  return(values)
}


#' Read 3d array (Analyze/Gate)
#' 
#' Read file with format "*.img, *.hdr" (Analyze/Gate) and return a \code{values} object.
#' 
#' @param file.name the file name
#' @param variable name of the variable to be assigned to the 3d data
#' @param voxel.origin coordinates of the "origin" voxel, i.e. the voxel cooresponding to indices (1,1,1)
#' 
#' @family R/W Arrays
#' @import oro.nifti
#' @export
read.3d.hdr <- function(file.name, variable='Dose[Gy]', voxel.origin=c(0,0,0))
{
  an <- readANALYZE(file.name)
  
  # array
  # check dimensioni 
  if(length(dim(an@.Data))==3) {Values.3d <- an@.Data}
  else if (length(dim(an@.Data))==4) {Values.3d <- an@.Data[,,,1]}
  else {stop('error in matrix dimension: ', length(dim(an@.Data)))}
  Nx <- an@dim_[2]
  Ny <- an@dim_[3]
  Nz <- an@dim_[4]
  Nv <- 1
  
  # coordinates
  dx <- an@pixdim[2]
  dy <- an@pixdim[3]
  dz <- an@pixdim[4]
  x <- seq(from=voxel.origin[1], by=dx, length.out=Nx)
  y <- seq(from=voxel.origin[2], by=dy, length.out=Ny)
  z <- seq(from=voxel.origin[3], by=dz, length.out=Nz)
  
  # crea oggetto values
  values <- list(values=Values.3d, x=x, y=y, z=z, Nx=Nx, Ny=Ny, Nz=Nz, Nv=Nv, variables=variable, file=file.name)
  class(values) <- 'values'
  return(values)
}


#' Read 3d DICOM 
#' 
#' Read 3D image from DICOM dile and return a \code{values} object.
#' Note: it assumes a regurlary spaced slices, constant pixel dimensions among slices
#' and the same ImagePositionPatient for each slice.
#' 
#' @param dicom.folder Path name to the DICOM directory.
#' @param recursive Search recursively down from the given path name.
#' @param exclude Exclude file names containing this character string.
#' @param verbose Flag to provide text-based progress bar.
#' @param variable name of the variable to be assigned to the 3d data.
#' @return A values object.
#' 
#' @family R/W Arrays, DICOM
#' @import oro.dicom
#' @export
read.3d.dicom <- function(dicom.folder, exclude=NULL, recursive=TRUE, verbose=TRUE, invert=TRUE, variable='HounsfieldNumber')
{
  message('New dicom implementation.')
  dcmImages <- readDICOM(path = dicom.folder, verbose = verbose,
                         recursive = recursive, exclude = exclude)
  dcm.info <- dicomTable(dcmImages$hdr)
  message("scaling images...")
  rs <- as.numeric(dcm.info$`0028-1053-RescaleSlope`)
  ri <- as.numeric(dcm.info$`0028-1052-RescaleIntercept`)
  if(length(rs)>0 & length(ri)>0) {
    for (i in 1:length(rs)) {
      dcmImages$img[[i]] <- dcmImages$img[[i]] * rs[i] + ri[i]
    }
  } else {
    message('slope & intercept not present in dicom, using dose grid scaling...')
    dose.grid.scaling <- as.numeric(dcm.info$`3004-000E-DoseGridScaling`)
    dcmImages$img[[1]] <- dcmImages$img[[1]] * dose.grid.scaling
  }

  if (invert & length(dcmImages$img)>1) {
    message("inverting images...")
    for (i in 1:length(dcmImages$img)) {
      dcmImages$img[[i]] <- dcmImages$img[[i]][seq(dim(dcmImages$img[[i]])[1], 1), ]
    }
  }

  if(length(dcmImages$img)>1) {
    message("sorting images...")
    zz <- as.numeric(dcm.info$`0020-1041-SliceLocation`)
    oz <- order(zz)
    imgs <- list()
    for (i in 1:length(oz)) {
      imgs[[i]] <- dcmImages$img[[oz[i]]]
    }
    dcmImages$img <- imgs
    Values.3d <- create3D(dcmImages)
  } else {
    Values.3d <- dcmImages$img[[1]]
    # inverte x - y
    Values.3d <- aperm(Values.3d, c(2, 1, 3))
    Values.3d <- DescTools::Rev(Values.3d, margin=2)
  }

  Nx <- dim(Values.3d)[1]
  Ny <- dim(Values.3d)[2]
  Nz <- dim(Values.3d)[3]
  Nv <- 1
  dxy <- as.numeric(unlist(strsplit(paste(dcm.info$`0028-0030-PixelSpacing`,
                                          collapse = " "), split = " ", fixed = TRUE)))
  vo <- as.numeric(unlist(strsplit(paste(dcm.info$`0020-0032-ImagePositionPatient`,
                                         collapse = " "), split = " ", fixed = TRUE)))
  if(length(dcmImages$img)>1) {
    voxel.origin <- c(mean(vo[seq(1, length(dxy) - 2, by = 3)]),
                      mean(vo[seq(2, length(dxy) - 1, by = 3)]),
                      mean(vo[seq(3, length(dxy), by = 3)]))
  } else {
    voxel.origin <- vo
  }
  message("origin: ", voxel.origin[1], ", ", voxel.origin[2],
          ", ", voxel.origin[3])
  dx <- mean(dxy[seq(1, length(dxy) - 1, by = 2)])
  dy <- mean(dxy[seq(2, length(dxy), by = 2)])
  x <- seq(from = voxel.origin[1], by = dx, length.out = Nx)
  y <- seq(from = voxel.origin[2], by = dy, length.out = Ny)
  if(length(dcmImages$img)>1) {
    z <- sort(zz) - min(zz) + voxel.origin[3]
    z <- sort(zz)
  } else {
    dz <- as.numeric(dcm.info$`0018-0050-SliceThickness`)
    z <- seq(from = voxel.origin[3], by = dz, length.out = Nz)
  }

  values <- list(values = Values.3d, x = x, y = y, z = z,
                 Nx = Nx, Ny = Ny, Nz = Nz, Nv = Nv, variables = variable)
  class(values) <- "values"
  message('created values with variable = ', values$variables)
  return(values)
}


#' Read 1D array
#' 
#' Read a file "*.1d" (PlanKIT format) and return a dataframe.
#' 
#' @param file.name the file name
#' @return a dataframe
#' 
#' @family R/W Arrays
#' @export
read.1d <- function(file.name) {
  
  #library(reshape)
  
  cat('reading values:', file.name, '\n')
  
  # apri connsessione
  file.1d <- file(file.name, "rb") # read binary
  
  # leggi l'header
  cat('reading header...\n')
  myline <- readLines(file.1d, n=4)
  
  # parsing dell'header
  # splitta la stringa in sottostringhe delimitate da uno o piu' spazi (regexpr: " +")
  myline.splitted <- unlist(strsplit(myline[1], ' +'))
  Nx <- as.numeric(myline.splitted[1])
  x.variable <- myline.splitted[2]
  x.type <- myline.splitted[3]
  myline.splitted <- unlist(strsplit(myline[2], ' +'))
  x <- as.numeric(myline.splitted)
  
  myline.splitted <- unlist(strsplit(myline[3], ' +'))
  Nv <- as.numeric(myline.splitted[1])
  variables <- myline.splitted[2:length(myline.splitted)]
  
  myline.splitted <- unlist(strsplit(myline[4], ' +'))
  values.format <- myline.splitted[1]
  
  Ntot <- Nx * Nv
  cat('number of elements:', Nv, 'x', Nx,  '=', Ntot, '\n')
  cat('independent variable:', x.variable, '\n')
  cat('dependent variables:', variables, '\n')
  
  # trasforma intervalli in coordinate puntuali
  if (x.type=='INTERVAL') {
    x.min <- x[1:Nx]
    x.max <- x[2:(Nx+1)]
    x <- (x[1:Nx] + x[2:(Nx+1)])/2
  }
  
  # crea e legge array
  if(values.format=='ASCII') {
    cat('Reading value data...\n')
    values <- t(array(as.numeric(unlist(strsplit(readLines(file.1d), ' +'))), dim=c(Nv, Nx)))
    #print(values)
    #array.1d <- array(readLines(file.1d), dim=c(Nv, Nx))
  }
  else if(values.format=='BINARY'){
    cat('warning: BINARY read still to be implemented.\n')
    return()
  }
  
  # chiudi file
  close(file.1d)
  
  # struttura data.frame
  values <- as.data.frame(values)
  names(values) <- variables
  values <- cbind(x, values)
  names(values)[1] <- x.variable
  if(x.type=='INTERVAL') {
    values <- cbind(x.min, x.max, values)
    names(values)[1:2] <- c(paste(x.variable, c('.min', '.max'), sep=''))
  }
  
  return(values)
}


#' Write 3D array (Pure-dek)
#' 
#' Write a \code{values} object in a file using the pure-dek format (*.3d)
#' 
#' @param values the \code{values} object
#' @param file.name the file name
#' 
#' @family R/W Array
#' @export
write.3d <- function(values, file.name) {
  
  Nx <- values[['Nx']]
  Ny <- values[['Ny']]
  Nz <- values[['Nz']]
  Nv <- values[['Nv']]
  variables <- values[['variables']]
  
  x <- values[['x']]
  y <- values[['y']]
  z <- values[['z']]
  
  dx <- (x[Nx]-x[1])/(Nx-1) # assume spaziatura costante
  dy <- (y[Ny]-y[1])/(Ny-1) # assume spaziatura costante
  dz <- (z[Nz]-z[1])/(Nz-1) # assume spaziatura costante
  
  x.int <- c(x-dx/2, x[Nx]+dx/2)
  y.int <- c(y-dy/2, y[Ny]+dy/2)
  z.int <- c(z-dz/2, z[Nz]+dz/2)
  
  
  # scrive header
  con <- file(file.name, "w") # open for writing in text mode
  writeLines(paste(Nx, "  X[mm]   INTERVAL\n"), con=con, sep='')
  writeLines(paste(x.int, collapse=' '), con=con)
  writeLines(paste(Ny, "  Y[mm]   INTERVAL\n"), con=con, sep='')
  writeLines(paste(y.int, collapse=' '), con=con)
  writeLines(paste(Nz, "  Z[mm]   INTERVAL\n"), con=con, sep='')
  writeLines(paste(z.int, collapse=' '), con=con)
  writeLines(paste(Nv, variables, collapse=' '), con=con)
  writeLines('BINARY', con=con)
  close(con)
  
  # scrive dati binari
  con <- file(file.name, "ab") # open for writing in text mode
  writeBin(as.double(values[['values']]),con=con)
  close(con)
}


# ANALIZE FORMAT ---------------------------------------------------------------

#' Write 3D array (analyze)
#' 
#' Write a \code{values} object in two files using the analyze format (*.hdr, *.img)
#' 
#' Note: with the analyze format only a single variable array from \code{values} can be saved in the same files.
#' 
#' @param values the \code{values} object
#' @param variable the variable name (default: first variable in values).
#' @param file.name the file name. Two files will be produced: \code{<file.name>.hdr} and \code{<file.name>.img}
#' @param gzipped compress the file (boolean, optional)
#' 
#' @family Analyze Format
#' @export
#' @import oro.nifti
write.analyze <- function(values, variable=NULL,
                          file.name, gzipped=FALSE)
{

  # librerie
  #library(oro.nifti)

  Nx <- values$Nx
  Ny <- values$Ny
  Nz <- values$Nz
  dx <- mean(values$x[2:Nx]-values$x[1:(Nx-1)]) # assume spaziatura costante
  dy <- mean(values$y[2:Ny]-values$y[1:(Ny-1)]) # assume spaziatura costante
  dz <- mean(values$z[2:Nz]-values$z[1:(Nz-1)]) # assume spaziatura costante

  # estrai 3d array
  if(is.null(variable)) {variable <- values$variables[1]}
  if(values$Nv > 1) {
    iv <- which(values$variable)
    array.3d <- values$values[iv,,,]
  } else {
    array.3d <- values$values
  }

  # crea oggetto analyze
  analyze.3d <- anlz(array.3d, datatype=16) # datatype 16 è "float" (http://andreaattili.org/appunti/?p=446) (gate accetta solo float o unsigned short)
  analyze.3d@"pixdim" <- c(0, dx, dy, dz, 0, 0, 0, 0)

  message('writing image...')
  print(analyze.3d)

  # salva oggetto analyze
  writeANALYZE(analyze.3d, file.name, gzipped=gzipped)
}


# UTILITIES --------------------------------------------------------------------


#' Sanitize file name
#' 
#' @param filename the file name
#' @return sanitized file name
#' 
#' @family Utilities
#' @export
sanitize.filename <- function(filename)
{
  filename <- gsub('%', '_', filename)
  return(filename)
}

# DICOM ------------------------------------------------------------------------

#' Convert DICOM header file/data.frame to nested list
#' 
#' Convert the DICOM header file or the header plain data frame obtained from the oro.dicom::readDICOMfile()
#' to a nested list of data. The nested list is analogous to the structure obtained in MatLab using
#' the function dicominfo() and it is easier to navigate and parse than the plain data frame.
#' 
#' Note: info about the element types
#' can be found at http://dicom.nema.org/dicom/2013/output/chtml/part05/sect_6.2.html
#' 
#' @param file the DICOM header file (alternative to the header data frame).
#' @param header the header data frame (alternative to DICOM header file).
#' @param index the starting row of the header data.frame to be parsed.
#' @param lev starting nesting level.
#' @param ... parameters passed to readDICOMFile().
#' 
#' @family DICOM
#' @export
#' @import oro.dicom
header2list <- function(file=NULL, header=NULL, index=1, lev=0, ...) {
  
  if(!is.null(file)) {
    message('reading DICOM file...')
    header <- readDICOMFile(fname = file, pixelData = FALSE, ...)[['hdr']]
  }
  
  if(index==1) {
    message('nesting DICOM data...')
  }
  
  my.list <- list()
  Nrow <- nrow(header)
  lev <- lev + 1
  item <- 1
  
  while(index <= Nrow) {
    
    # item
    if(header[['name']][index] == 'Item') {
      item.name <- paste0('Item_', item)
      item <- item + 1
      index <- index + 1
      
      temp <- header2list(file=NULL, header=header, index=index, lev=lev)
      my.list[[item.name]] <- temp[[1]]
      index <- temp[[2]]
      
      #return(list(my.list, index))
    }
    
    # sequence
    else if(header[['value']][index] == 'Sequence') {
      sequence.name <- header[['name']][index]
      index <- index + 1
      
      temp <- header2list(file=NULL, header=header, index=index, lev=lev)
      my.list[[sequence.name]] <- temp[[1]]
      index <- temp[[2]]
      
      #return(list(my.list, index))
    }
    
    # esci 
    else if( (header[['name']][index] == 'SequenceDelimitationItem') | (header[['name']][index] == 'ItemDelimitationItem')) {
      index <- index + 1
      return(list(my.list, index))
    }
    
    # campo semplice
    else {
      item.name <- header[['name']][index]
      value <- header[['value']][index]
      #message(index, '-->', lev )
      #message(item.name, ' = ', value)
      my.list[[item.name]] <- value
      index <- index + 1
    }
  }
  if(lev>1) return(list(my.list, index)) else return(my.list)
}
