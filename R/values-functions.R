# VALUES -----------------------------------------------------------------------

#' Create values object
#' 
#' Create values object from array
#' 
#' @param array the 3D array
#' @param variable the name of the variable
#' @param x,y,z coordinate vectors
#' @return a \code{values} object
#' 
#' @family Values
#' @export
#' 
create.values <- function(array.values=NULL, variables='variable name', x, y, z)
{
  
  Nx <- length(x)
  Ny <- length(y)
  Nz <- length(z)
  Nv <- length(variables)
  
  if(is.null(array.values)) {
    if(Nv>1) {
      array.values <- array(0, dim=c(Nx, Ny, Nz, Nv))
    } else {
      array.values <- array(0, dim=c(Nx, Ny, Nz))
    }
  }
  
  values <- list()
  values$values <- array.values
  values$x <- x
  values$y <- y
  values$z <- z
  values$Nx <- Nx
  values$Ny <- Ny
  values$Nz <- Nz
  values$Nv <- Nv
  values$variables <- variables
  values$file <- NULL
  
  # check dimensioni
  my.dims <- dim(values$values)
  if(my.dims[1]!=values$Nx | my.dims[2]!=values$Ny | my.dims[3]!=values$Nz) {
    message('create.values warning: inconsistent size with coordinates')
  }
  
  return(values)
}


#' Get values object (PlanKIT)
#' 
#' Get the values matrices for a plan (PlanKIT).
#' @param plan The plan object.
#' @family Values
#' @export
get.values <- function(plan)
{
  if(is.null(plan[['outputValuesFile']])) {
    cat('The plan "', plan[['name']], '" has no values file.\n', sep='')
    return(NULL)
  } 
  return(read.3d( paste(plan[['name']], '/', plan[['outputValuesFile']], sep='') ))
}


#' Get values object (Gate)
#' 
#' Get the values matrices for a plan (PlanKIT).
#' @param plan.gate The plan object.
#' @family Values
#' @export
get.values.gate <- function(plan.gate)
{
  values <- list()
  variables <- plan.gate$computingValues.gate
  for(i in 1:length(variables))
  {
    variable.file <- file.variable.gate(variables[i])
    variable.file <- paste(plan.gate$name, '/output/', variable.file, sep='')
    values[[i]] <- read.3d.hdr(file.name=variable.file, variable=variables[i])
  }
  return(merge.values(values.list=values))
}


#' recupera matrice di valori in formato dataframe
#' 
#' @family Values
#' @export
get.values.dataframe <- function(plan)
{
  if(is.null(plan[['outputValuesFile']])) {
    cat('The plan "', plan[['name']], '" has no values file.\n', sep='')
    return(NULL)
  } 
  return(read.3d( paste(plan[['name']], '/', plan[['outputValuesFile']], sep='') ))
}


#' genera struttura dataframe da oggetto values
#' 
#' l'estrazione viene fatta sui VOI specificati
#' ogni riga corrisponde ad un voxel.
#' nel dataframe ogni voxel è associato a un VOI.
#' voxel associati a più voi vengono ripetuti su più righe.
#' 
#' @family Values
#' @export
get.dataframe.from.values <- function(values=values, vois=vois, variables=NULL, rois=NULL)
{
  # identifica i vois da mettere nel dataframe
  if(is.null(rois)) {
    rois <- vois$vois
    index.r <- 1:length(rois)
  } else {
    index.r <- which(vois$vois %in% rois) # gli indici dei voi si riferiscono al vettore vois$vois
  }
  
  # identifica le variabili da mettere nel dataframe
  if(is.null(variables)) {
    variables <- values$variables
    index.v <- 1:values$Nv
  } else {
    index.v <- which(values$variables %in% variables) # gli indici delle variabili si riferiscono al vettore values$variables
  }
  
  # crea le matrici X,Y,Z
  X <- array(values$x, dim=c(values$Nx,values$Ny,values$Nz))
  Y <- array(values$y, dim=c(values$Ny,values$Nx,values$Nz))
  Y <- aperm(Y, c(2,1,3))
  Z <- array(values$z, dim=c(values$Nz,values$Nx,values$Ny))
  Z <- aperm(Z, c(2,3,1))
  
  # crea dataframe per ogni voi specificato nella lista rois
  for(i in 1:length(index.r))
  {
    message('extracting dataframe for: ', vois$vois[index.r[i]])
    
    index.voi <- get.voi.logical(vois, voi=vois$vois[index.r[i]])
    
    # estrai coordinate dei voi
    x <- X[index.voi]
    y <- Y[index.voi]
    z <- Z[index.voi]
    
    # ciclo sulle variabili
    for(j in 1:length(index.v)) {
      message('- variable: ', values$variables[index.v[j]])
      if(values$Nv > 1) {
        my.array <- values$values[index.v[j],,,][index.voi]
      } else {
        my.array <- values$values[,,][index.voi]
      }
      df.tmp <- data.frame(voi=vois$vois[index.r[i]], variable=values$variables[index.v[j]], x=x, y=y, z=z, value=my.array)
      if(i==1 & j==1) {df.all <- df.tmp} else {df.all <- rbind(df.all, df.tmp)}
    }
  }
  
  return(df.all)
}


#' recupera matrice sparsa da values
#' 
#' accetta values con più di una variabile. la selezione è fatta sulla prima.
#' (che può coincidere o meno con la variabile specificata)
#' 
#' @family Values
#' @export
get.sparse.array.from.values <- function(values, variable='Dose[Gy]', threshold=0)
{
  # identifica variabile
  if(values$Nv>1) {
    v.index <- which(values$variables==variable)
    my.array <- values$values[v.index,,,]
    voxel.index <- which(values$values[1,,,]>threshold)
  } else {
    my.array <- values$values
    voxel.index <- which(my.array>threshold)
  }
  
  # matrice sparsa
  sparse.array <- my.array[voxel.index]
  
  return(data.frame(voxel.index=voxel.index, value=sparse.array))
}


#' ricostruisce oggetto values 3d dalla matrice sparsa
#' 
#' l'oggetto values conterrà solo una variabile, con il nome specificato
#' 
#' @family Values
#' @export
#' 
get.values.from.sparse.array <- function(sparse.array, variable='Dose[Gy]', x, y, z)
{
  values <- list()
  
  Nx <- length(x)
  Ny <- length(y)
  Nz <- length(z)
  Nv <- 1
  
  values$values <- array(0, dim=c(Nx, Ny, Nz))
  values$values[sparse.array$voxel.index] <- sparse.array[,2]
  values$x <- x
  values$y <- y
  values$z <- z
  values$Nx <- Nx
  values$Ny <- Ny
  values$Nz <- Nz
  values$Nv <- Nv
  values$variables <- variable
  values$file <- ''
  
  return(values)
}


#' rimuove valori nulli non "fisici"
#' 
#' @family Values
#' @export
sanitize.values <- function(values)
{
  variables <- c('Alpha[Gy^(-1)]', 'Beta[Gy^(-2)]', 'RBE', 'DoseAveragedLET[keV/um]', 'DoseAveragedLET[keV/um]')
  for(vs in variables) {
    v <- which(values$variables==vs)
    if (length(v)>0) {
      vv <- values$values[v,,,]
      vv[vv<=0] <- NA
      values$values[v,,,] <- vv
    }
  }
  return(values)
}

#' Merge a list of values
#' 
#' Merge a list of values objects in a single object. It assumes that the variable names are unique.
#' @param values.list A list of values objects.
#' @return A values object
#' @family Values
#' @export
merge.values <- function(values.list)
{
  NV <- length(values.list)
  if(NV==1) {return(values.list[[1]])}
  values1 <- values.list[[1]]
  for(i in 2:NV) {
    values <- values.list[[i]]
    nv <- values$Nv
    for(v in 1:nv) {
      values1 <- add.array.values(values=values1, new.array=values$values, variable=values$variable[v])
    }
  }
  return(values1)
}


# ARRAY VALUES -----------------------------------------------------------------

#' Get origin coordinates
#' 
#' Get the coordinates of origin corner (corner of voxel [1,1,1]).
#' @param values the values object
#' @return a vector of the coordinates x, y and z.
#' @family Array Values
#' @export
get.origin <- function(values)
{
  x0 <- values$x[1] - (values$x[2] - values$x[1])/2
  y0 <- values$y[1] - (values$y[2] - values$y[1])/2
  z0 <- values$z[1] - (values$z[2] - values$z[1])/2
  return(c(x0, y0, z0))
}

#' Set origin coordinates
#' 
#' Set the coordinates of origin corner (corner of voxel [1,1,1]).
#' @param coordinates a vector of the coordinates x, y and z.
#' @return the translated values object
#' @family Array Values
#' @export
set.origin <- function(values, coordinates)
{
  x0 <- values$x[1] - (values$x[2] - values$x[1])/2
  y0 <- values$y[1] - (values$y[2] - values$y[1])/2
  z0 <- values$z[1] - (values$z[2] - values$z[1])/2
  
  values$x <- values$x - x0 + coordinates[1]
  values$y <- values$y - y0 + coordinates[2]
  values$z <- values$z - z0 + coordinates[3]
  
  return(values)
}

#' estrai array della variabile specificata dai values
#' 
#' @family Array Values
#' @export
get.array.values <- function(values, variable)
{
  index.v <- which(values$variables==variable) # identifica indice della variabile specificata
  if(length(index.v)==0) {message('Variable not found.'); return()}
  if(values$Nv > 1) {
    return(values$values[index.v,,,])
  } else {
    return(values$values)
  }
}


#' aggiungi array di una nuova variabile all'oggetto values
#' 
#' @family Array Values
#' @export
add.array.values <- function(values, new.array, variable='New Data')
{
  # check per vedere se la variabile specificata già esiste
  if(length(which(values$variables==variable)) >= 1) {
    message('Variable already in values')
    return()
  }
  
  if(values$Nv > 1) {
    v <- aperm(values$values, c(2,3,4,1))
  } else {
    v <- values$values
  }
  v <- array(c(v, new.array), c(values$Nx, values$Ny, values$Nz, values$Nv+1))
  values$values <- aperm(v, c(4,1,2,3))
  values$Nv <- values$Nv+1
  values$variables[values$Nv] <- variable
  
  return(values)
}


#' identifica slice corrispondente alla coordinata
#' 
#' @family Profiles
#' @export
identify.slice <- function(x, xvec)
{
  return(which( abs(x-xvec) == min(abs(x-xvec)))[1])
}


#' Get profile
#' 
#' Get a profile (values along a specific axis, or voxel line) from a \code{values} object. The axis is specifed trough the definition of two of the following three coordinates: x, y, z.
#' 
#' @param values the \code{values object}
#' @param variable the variable to be profiled
#' @param x,y,z definition of the axis (only two coordinates should be defined). For example \code{x=4.4} and \code{z=0.5} specifies a profile along the y-axis passing through x=4.4 and z=0.5
#' @param integrate perform an integration over the two coordinates (boolean, optional)
#' @return a dataframe containing the specified profile
#' 
#' @family Values
#' @export
get.profile <- function(values, variable=NULL, x=NULL, y=NULL, z=NULL, integrate=FALSE)
{
  
  # variable
  v <- NULL
  if(!is.null(variable)) {v <- which(values$variables==variable)}
  if(length(v)==0) {v <- 1}
  variable <- values$variable[v]
  message('using ', variable, ' ...')
  if(values$Nv>1) {values$values <- values$values[v,,,]}
  
  Nx <- Ny <- Nz <- NULL
  if(!is.null(x)) {Nx <- identify.slice(x, values$x)}
  if(!is.null(y)) {Ny <- identify.slice(y, values$y)}
  if(!is.null(z)) {Nz <- identify.slice(z, values$z)}
  
  # recupera coordinate
  if(is.null(x)) {coord.name <- 'x[mm]'; coord.value <- values$x}
  if(is.null(y)) {coord.name <- 'y[mm]'; coord.value <- values$y}
  if(is.null(z)) {coord.name <- 'z[mm]'; coord.value <- values$z}
  
  # calcoli preliminari per l'integrazione
  if(integrate) {
    variable <- paste(variable, '[mm^2]', sep='')
    
    dx <- diff(points2intervals(values$x))
    dy <- diff(points2intervals(values$y))
    dz <- diff(points2intervals(values$z))

  }
  
  # calcola profili
  if(!is.null(Nx) & !is.null(Ny)) {
    if(integrate) {
      profile.value <- apply(values$values, c(3), sum) * dx * dy
    } else {
      profile.value <- values$values[Nx,Ny,]
    }
  }
  if(!is.null(Nx) & !is.null(Nz)) {
    if(integrate) {
      profile.value <- apply(values$values, c(2), sum) * dx * dz
    } else {
      profile.value <- values$values[Nx,,Nz]
    }
  }
  if(!is.null(Ny) & !is.null(Nz)) {
    if(integrate) {
      profile.value <- apply(values$values, c(1), sum) * dy * dz
    } else {
      profile.value <- values$values[,Ny,Nz]
    }
  }
  

  
  profile.df <- data.frame(variable=variable, axis=coord.name, depth=coord.value, value=profile.value)
  
  return(profile.df)
}



# UTILITIES --------------------------------------------------------------------

#' Evaluate intervals
#' 
#' Evaluate intervals vector from points vector.
#' 
#' Voxel coordinates are usually specified by indicating (x,y,z) the middle point. Another way is to indicate the intervals (xmin,xmax), (ymin,ymax) (zmin,zmax). Since it is assumed that voxels are contiguous, xmax[i] = xmin[i-1], for each axis the coordinates can be stored in a single vector.
#' 
#' @param x the puntual coordinate vector
#' @return the intervals vector (the length of this vector is Nx+1, where Nx = length(x))
#' @export
#' @family Utilities
points2intervals <- function(x)
{
  if(length(x)>1) {
    dx <- diff(x)
    xi <- c(x[1]-mean(dx)/2,
            x[1:(length(x)-1)]+dx/2,
            x[length(x)]+mean(dx)/2)
  } else {xi <- c(x-0.5, x+0.5)}
  return(xi)
}


#' Removes values outside VOI
#' 
#' @param values The original values
#' @param vois The vois object
#' @param voi the name of the VOI
#' @param index.voi It is possible to pass directly the logical index array.
#' @return A values object
#' @export
#' @family Vois, Values, Utilities
remove.values.outside.voi <- function(values, vois, voi, index.voi=NULL)
{
  
  if(is.null(index.voi)) {
    index.voi <- get.voi.logical(vois, voi)
  }
  
  Nv <- values$Nv
  if(Nv==1) {
    values$values[!index.voi] <- NA
  } else {
    for(i in 1:Nv) {
      values.tmp <- values$values[i,,,]
      values.tmp[!index.voi] <- NA
      values$values[i,,,] <- values.tmp
    }
  }
  
  return(values)
}

#' Create values form events
#' 
#' Uses the (x,y,z) coordinates of a shower of events to create a values object with the 3D histogram of the events.
#' 
#' @param Xe event x coorditate vector.
#' @param Ye event y coordinate vector.
#' @param Ze event z coordinate vector.
#' @param weight weight vector.
#' @param xbin bins x coordinates (extremes).
#' @param ybin bins y coordinates (extremes).
#' @param zbin bins z coordinates (extremes).
#' @param variable The variable name.
#' @return A values object
#' @exports
#' @family Values, Utilities
generate.values.from.events <- function(Xe, Ye, Ze, xbin, ybin, zbin, weight=1, variable='Activity', group=10000)
{
  
  message('histogramming the ', variable, '...')
  if(length(weight)==1) {weight=rep(weight, length(Xe))}
  
  # intervalli (da coordinate puntuali)
  Nx <- length(xbin) - 1
  Ny <- length(ybin) - 1
  Nz <- length(zbin) - 1
  x <- (xbin[1:Nx] + xbin[2:(Nx+1)])/2
  y <- (ybin[1:Ny] + ybin[2:(Ny+1)])/2
  z <- (zbin[1:Nz] + zbin[2:(Nz+1)])/2
  
  # cuts (se riuscissi a fare direttamente un cut=voxel index...)
  Xe.c <- cut(Xe, breaks=xbin)
  Ye.c <- cut(Ye, breaks=ybin)
  Ze.c <- cut(Ze, breaks=zbin)
  voxel.index <- as.numeric(interaction(Xe.c, Ye.c, Ze.c)) # questo è quasi un miracolo della fede.
  
  # crea oggetto values (vuoto)
  # values <- create.values(variables=variable, x=x, y=y, z=z)
  
  #XYZ <- aggregate(list(weigth=weight), list(x=Xe.c, y=Ye.c, z=Ze.c), sum)
  XYZ <- aggregate(list(value=weight), list(voxel.index=voxel.index), sum)
  
  return(get.values.from.sparse.array(sparse.array=XYZ, variable=variable, x=x, y=y, z=z))
}
  
