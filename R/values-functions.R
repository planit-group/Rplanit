# VALUES -----------------------------------------------------------------------

#' Create values object
#' 
#' Create values object from array. Note: the value distributions are stored in a multidimensional array:
#' values$values with dim(values$values) = c(Nv, Nx, Nz, Ny).
#' 
#' @param array the 3D array. Optional. If it is not specified an array of zeros is assumed.
#' @param variables The name of the variable (can be a vector of names).
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
      array.values <- array(0, dim=c(Nv, Nx, Ny, Nz))
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
  if(length(my.dims)==3) {
    if(my.dims[1]!=values$Nx | my.dims[2]!=values$Ny | my.dims[3]!=values$Nz) {
      message('create.values warning: inconsistent size with coordinates')
    }
  } else {
    if(my.dims[1]!=values$Nv | my.dims[2]!=values$Nx | my.dims[3]!=values$Ny | my.dims[4]!=values$Nz) {
      message('create.values warning: inconsistent size with coordinates')
    }
  }
  
  return(values)
}


#' Get values object
#' 
#' Get the values matrices for a plan.
#' @param plan The plan object.
#' @family Values
#' @export
get.values <- function(plan, ...) UseMethod("get.values")

#' Get values object (PlanKIT)
#' 
#' Get the values matrices for a plan (PlanKIT).
#' @param plan The plan object.
#' @family Values
#' @export
get.values.plankit.plan <- function(plan)
{
  if(is.null(plan[['outputValuesFile']])) {
    cat('The plan "', plan[['name']], '" has no values file.\n', sep='')
    return(NULL)
  } 
  return(read.3d( paste(plan[['name']], '/', plan[['outputValuesFile']], sep='') ))
}

#' Get values object (Gate)
#' 
#' Get the values matrices for a plan (Gate).
#' @param plan.gate The plan object.
#' @param center.coordinates The coordinates of the center of the values 3D distribution. If it is NULL, the center of the CT is assumed.
#' @family Values
#' @export
get.values.gate.plan <- function(plan, center.coordinates=NULL)
{
  values <- list()
  variables <- plan$computingValues
  for(i in 1:length(variables))
  {
    variable.file <- file.variable.gate(variables[i])
    variable.file <- paste(plan$name, '/output/', variable.file, sep='')
    values[[i]] <- read.3d.hdr(file.name=variable.file, variable=variables[i])
  }
  values <- merge.values(values.list=values)
  
  # recupera le coordinate assolute dalla CT
  if(is.null(center.coordinates)) {
    ct <- get.ct(plan)
    center.coordinates[1] <- sum(range(ct$x))/2
    center.coordinates[2] <- sum(range(ct$y))/2
    center.coordinates[3] <- sum(range(ct$z))/2
    print(center.coordinates)
  }
  
  # traslazione del centro
  x.c <- sum(range(values$x))/2
  y.c <- sum(range(values$y))/2
  z.c <- sum(range(values$z))/2
  values$x <- values$x + (center.coordinates[1]-x.c)
  values$y <- values$y + (center.coordinates[2]-y.c)
  values$z <- values$z + (center.coordinates[3]-z.c)
  
  return(values)
}


#' Get sparse array values (Gate)
#' 
#' Get the sparse array evaluated for a plan (Gate).
#' Note: the plan has to be evaluated via run.gate.forward() with the option evaluate.sparse.arrays=TRUE.
#' @param plan.gate The plan object.
#' @family Values
#' @export
get.sparse.array.gate <- function(plan.gate)
{
  load(paste(plan.gate$name, 'values.sparse.Rdata', sep='/'))
  return(values.sparse)
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
    
    for(i in 1:values$Nv) {
      my.array <- values$values[i,,,]
      sparse.array <- my.array[voxel.index]
      if(i==1) {
        df <- data.frame(voxel.index=voxel.index, value=sparse.array, variable=as.factor(values$variables[i]))
      } else {
        df <- rbind(df, data.frame(voxel.index=voxel.index, value=sparse.array, variable=as.factor(values$variables[i])))
      }
    }
    return(df)
    
  } else {
    my.array <- values$values
    voxel.index <- which(my.array>threshold)
    sparse.array <- my.array[voxel.index]
    return(data.frame(voxel.index=voxel.index, value=sparse.array))
  }
  
}


#' Create values from sparse array
#' 
#' Create a values object from a sparse array.
#'
#' @param sparse.array The sparse array. Note: the sparse array is a dataframe containing at least the following columns: 'voxel.index' and 'value'. Optionally it could contain also the column 'variable'. The latter is mandatory if there are multiple variables.
#' @param variables The name of the variable(s). If there is only one variable specified and the sparse array dose not contain the column "variable", it is assumed that all the value correspond to the specified variable. It it is NULL, the values object will contains all the variable found in the sparse array.
#' @return Values object.
#' Note: the sparse array is a dataframe containing at least the following columns: 'voxel.index' and 'value'
#' @family Values
#' @export
#' 
get.values.from.sparse.array <- function(sparse.array, variables=NULL, x, y, z)
{
  values <- list()
  
  Nx <- length(x)
  Ny <- length(y)
  Nz <- length(z)
  
  # variables
  if(is.null(variables)) {
    if('variable' %in% colnames(sparse.array)) {
      variables <- as.character(unique(sparse.array$variable))
    } else {
      warning('variable not specified, using dummy...')
      variables <- 'dummy'
    }
  } else {
    if('variable' %in% colnames(sparse.array)) {
      sparse.array <- subset(sparse.array, variable %in% variables)
      variables <- as.character(sparse.array$variable)
    } else {
      if(length(variables)>1) {stop('variable names not present in sparse array')}
    }  
  }
  Nv <- length(variables)
  
  values <- create.values(variables=variables, x=x, y=y, z=z)
  
  if(Nv==1) {values$values[sparse.array$voxel.index] <- sparse.array[,2]}
  else {
    for(i in 1:Nv) {
      temp.array <- array(0, dim=c(Nx, Ny, Nz))
      temp.array[sparse.array$voxel.index] <- sparse.array$value[sparse.array$variable==variables[i]]
      values$values[i,,,] <- temp.array
    }
  }
  
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


# VALUES MANIPULATION ----------------------------------------------------

#' Get origin coordinates
#' 
#' Get the coordinates of origin corner (corner of voxel [1,1,1]).
#' @param values the values object
#' @return a vector of the coordinates x, y and z.
#' @family Values Manipulation
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
#' @family Values Manipulation
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
#' @family Values Manipulation
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
#' @family Values Manipulation
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


# PROFILES ---------------------------------------------------------------------

#' Ray-tracing (intersections)
#' 
#' Evaluate the intersections along voxels 
#' @export
#' @family Values, Profiles
ray.tracing <- function(ray, xi, yi, zi, values=NULL) {
  
  if(!is.null(values)) {
    # coordinate piani tra i voxels
    xi <- create.intervals(values$x)
    yi <- create.intervals(values$y)
    zi <- create.intervals(values$z)   
  }
  
  R <- c(ray$X, ray$Y, ray$Z)
  rn <- c(ray$xn, ray$yn, ray$zn)
  px <- py <- pz <- list()
  
  # calcola intersezione dei piani
  if(rn[1]!=0) {
    tx <- (xi - R[1])/rn[1]
    px.x <- R[1]+tx*rn[1]
    px.y <- R[2]+tx*rn[2]
    px.z <- R[3]+tx*rn[3]
    px <- data.frame(x=px.x, y=px.y, z=px.z, t=tx)
  } else {
    tx <-NULL
    px <- data.frame(x=NA, y=NA, z=NA, t=NA)
  }
  if(rn[2]!=0) {
    ty <- (yi - R[2])/rn[2]
    py.x <- R[1]+ty*rn[1]
    py.y <- R[2]+ty*rn[2]
    py.z <- R[3]+ty*rn[3]
    py <- data.frame(x=py.x, y=py.y, z=py.z, t=ty)  
  } else {
    ty <- NULL
    py <- data.frame(x=NA, y=NA, z=NA, t=NA)
  }
  if(rn[3]!=0) {
    tz <- (zi - R[3])/rn[3]
    pz.x <- R[1]+tz*rn[1]
    pz.y <- R[2]+tz*rn[2]
    pz.z <- R[3]+tz*rn[3]
    pz <- data.frame(x=pz.x, y=pz.y, z=pz.z, t=tz)    
  } else {
    tz <- NULL
    pz <- data.frame(x=NA, y=NA, z=NA, t=NA)
  }
  
  # merge di punti
  p <- rbind(px, py, pz)
  
  # elimina punti fuori array
  p <- subset(p, x>=min(xi, na.rm=TRUE) & x<=max(xi, na.rm=TRUE) &
                y>=min(yi, na.rm=TRUE) & y<=max(yi, na.rm=TRUE) &
                z>=min(zi, na.rm=TRUE) & z<=max(zi, na.rm=TRUE))
  
  # sorting dei punti
  p <- p[order(p$t),]
  
  return(p)
}


#' Get profile
#' 
#' Get a profile (values along a specific axis, or voxel line) from a \code{values} object. The axis is specifed trough the definition of two of the following three coordinates: x, y, z, or with a single ray object for arbitrary directions.
#' 
#' @param values the \code{values object}
#' @param variable the variable to be profiled. If it is not specified, the fist available variable will be used.
#' @param x,y,z definition of the axis (only two coordinates should be defined). For example \code{x=4.4} and \code{z=0.5} specifies a profile along the y-axis passing through x=4.4 and z=0.5
#' @param integrate perform an integration over the two coordinates (boolean, optional). Note: integration is not yet implemented around an arbitrary ray.
#' @param ray The ray object, Note the ray is a data.frame defined by 6 components (X, Y, Z, xn, yn, zn), i.e. point coordinates + normalized vector direction.
#' @param return.voxel.index Returns also the index of the voxel crossed by the ray (used only for a ray profile).
#' @param return.xyz Returns also the x,y,z of the points along the ray (used only for a ray profile).
#' @return a dataframe containing the specified profile.
#' @family Values, Values Manipulation, Profiles
#' @todo Generation of multi-variable profiles (a list?). If no variable is specified all the available variable will be used. Integration around arbitrary ray. Gestisci in maniera uniforme i casi ortogonali rispetto al caso generale "ray".
#' @export
get.profile <- function(values, variable=NULL, x=NULL, y=NULL, z=NULL, integrate=FALSE, ray=NULL, return.voxel.index=FALSE, return.xyz=FALSE)
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
  
  # profilo su direzione arbitraria (ray-tracing)
  if(!is.null(ray)) {
    
    # coordinate piani tra i voxels
    #dx <- mean(diff(values$x)); dy <- mean(diff(values$y)); dz <- mean(diff(values$z))
    xi <- create.intervals(values$x)
    yi <- create.intervals(values$y)
    zi <- create.intervals(values$z)
    
    # evaluate ray-tracing
    p <- ray.tracing(ray, xi, yi, zi)
    
    np <- nrow(p)-1
    
    # punti centrali dei segmenti
    p.mid <- data.frame(x=(p$x[1:np]+diff(p$x)/2), y=(p$y[1:np]+diff(p$y)/2), z=(p$z[1:np]+diff(p$z)/2), t=(p$t[1:np]+diff(p$t)))
    
    # identifica voxels
    voxel.index <- get.voxel.index(Xv=values$x, Yv=values$y, Zv=values$z, x=p.mid$x, y=p.mid$y, z=p.mid$z)
    
    coord.name <- 'r [mm]'
    coord.value <- p.mid$t - min(p$t) # la coordinata è zero all'ingresso
    profile.value <- values$values[voxel.index]
  }

  
  profile.df <- data.frame(variable=variable, axis=coord.name, depth=coord.value, value=profile.value, variable='variable')
  
  if(return.voxel.index & !is.null(ray)) {profile.df$voxel.index <- voxel.index}
  if(return.xyz  & !is.null(ray)) {
    profile.df$x <- p.mid$x
    profile.df$y <- p.mid$y
    profile.df$z <- p.mid$z
  }
  
  return(profile.df)
}



# VALUES UTILITIES -------------------------------------------------------------

#' Tri-linear interpolation
#' 
#' @param values Values object
#' @param x,y,z cartesian coordinate on which evaluate the interpolation.
#' @return a values object (with elements at coordinates x,y,z).
#' 
#' @family ValuesUtilities
#' @export
interpolate.values <- function(values, x, y, z)
{
  Nx <- length(x); Ny <- length(y); Nz <- length(z)
  
  # crea oggetto values
  values.intrp <- create.values(variables=values$variables, x=x, y=y, z=z)
  
  for(v in 1:values$Nv) {
    if(values$Nv==1) {array.v <- values$values} else {array.v <- values$values[v,,,]}
    array.x <- array(0, dim=c(Nx, values$Ny, values$Nz))
    array.y <- array(0, dim=c(Nx, Ny, values$Nz))
    array.z <- array(0, dim=c(Nx, Ny, Nz))
    
    message('interpolating ', values$variables[v])
    
    # interpolazione lungo x
    for(iy in 1:values$Ny) {
      for(iz in 1:values$Nz) {
        vx.fun <- approxfun(x=values$x, y=array.v[,iy,iz], rule=2)
        vx <- vx.fun(x)
        array.x[,iy,iz] <- vx
      }
    }
    
    # interpolazione lungo y
    for(ix in 1:Nx) {
      for(iz in 1:values$Nz) {
        vy.fun <- approxfun(x=values$y, y=array.x[ix,,iz], rule=2)
        vy <- vy.fun(y)
        array.y[ix,,iz] <- vy
      }
    }
    
    # interpolazione lungo z
    for(ix in 1:Nx) {
      for(iy in 1:Ny) {
        vz.fun <- approxfun(x=values$z, y=array.y[ix,iy,], rule=2)
        vz <- vz.fun(z)
        array.z[ix,iy,] <- vz
      }
    }
   
    
    if(values$Nv==1) {
      values.intrp$values <- array.z
    } else {
      values.intrp$values[v,,,] <- array.z
    }
  }
  
  return(values.intrp)
}


#' identifica slice corrispondente alla coordinata
#' 
#' @family ValuesUtilities
#' @export
identify.slice <- function(x, xvec)
{
  return(which( abs(x-xvec) == min(abs(x-xvec)))[1])
}

#' Evaluate intervals
#' 
#' Evaluate intervals vector from points vector.
#' 
#' Voxel coordinates are usually specified by indicating (x,y,z) the middle point. Another way is to indicate the intervals (xmin,xmax), (ymin,ymax) (zmin,zmax). Since it is assumed that voxels are contiguous, xmax[i] = xmin[i-1], for each axis the coordinates can be stored in a single vector.
#' 
#' @param x the puntual coordinate vector
#' @return the intervals vector (the length of this vector is Nx+1, where Nx = length(x))
#' @export
#' @family Values Utilities
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
#' @family Vois, Values, Values Utilities
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
#' @param sparse.matrix Return a sparse matrix if TRUE.
#' @return Values object or a sparse matrix.
#' @export
#' @family Values, Values Utilities
generate.values.from.events <- function(Xe, Ye, Ze,
                                        xbin, ybin, zbin,
                                        weight=1, variable='Activity', sparse.matrix=FALSE)
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

  if(sparse.matrix) {
    return(XYZ)
  } else {
    return(get.values.from.sparse.array(sparse.array=XYZ, variables=variable, x=x, y=y, z=z))
  }
}
  
#' Generate bin intervals from points
#' 
#' @param x A vector.
#' @return A vector (with N+1 elements)
#' @family Values, Values Utilities
#' @export
create.intervals <- function(x)
{
  dx <- diff(x)
  nx <- length(x)
  xi <- c(x[1] - dx[1]/2, x[1:(nx-1)]+dx/2, x[nx]+dx[nx-1]/2)
  return(xi)
}

#' Get voxel index from coordinates
#' 
#' @param x,y,z coordinate vectors.
#' @param Xv,Yv,Zv coordinates of the voxels.
#' @return A vector.
#' @family Values, Values Utilities
#' @export
get.voxel.index <- function(x, y, z, Xv, Yv, Zv)
{
  xbin <- create.intervals(Xv)
  ybin <- create.intervals(Yv)
  zbin <- create.intervals(Zv)
  
  X.c <- cut(x, breaks=xbin)
  Y.c <- cut(y, breaks=ybin)
  Z.c <- cut(z, breaks=zbin)
  return(as.numeric(interaction(X.c, Y.c, Z.c)))
}

#' Get coordinates from voxel index
#' 
#' @param x,y,z coordinate vectors corresponding to the 3D array.
#' @param voxel.index The voxels index (can be a vectors of voxel index).
#' @return A data.frame (x,y,z) for more than a single voxel index. 
#' @family Values, Values Utilities
#' @export
get.coordinates <- function(voxel.index, x, y, z)
{
  # crea le matrici X,Y,Z
  Nx <- length(x); Ny <- length(y); Nz <- length(z)
  # X
  A <- array(x, dim=c(Nx,Ny,Nz))
  X <- A[voxel.index]
  # Y
  A <- array(y, dim=c(Ny,Nx,Nz)); A <- aperm(A, c(2,1,3))
  Y <- A[voxel.index]
  # Z
  A <- array(z, dim=c(Nz,Nx,Ny)); A <- aperm(A, c(2,3,1))
  Z <- A[voxel.index]

  return(data.frame(x=X, y=Y, z=Z))
}
