# DVH --------------------------------------------------------------------------


#' Evaluate DVH
#' 
#' Evaluate the DVH for a single VOI and a specified variable.
#' 
#' @param values the \code{values} object
#' @param vois the \code{vois} object
#' @param voi the VOI name
#' @param variable the name of the variable (optional if \code{values} contains only one variable)
#' @param with.na leave explicitly NA values in DVH representing voxels were the specified variable is not defined (optional, boolean)
#' @param with0 put explicitly in the DVH a point at (0, max(volume)). It is possible to specify directly the min value to use for this DVH point, in this case if with.0 is a number !=0, the point will be at (with.0, max(volume)) (optional)
#' @param name Overwrite the default name derived from the voi name.
#' @return The \code{dvh} object; a list consisting of:
#' \item{value}{vector containing the value data}
#' \item{valume}{vector containing the volume fractions}
#' \item{volume.tot}{total volume of the VOI in mm^3}
#' \item{voi}{the name of the VOI}
#' \item{Nvoxel}{the number of voxels of the VOI}
#' 
#' @family DVH
#' @export
dvh.evaluate <- function(values=NULL, vois=NULL, voi='PTV', variable=NULL, with.na=FALSE, with.0=FALSE, name=NULL)
{
  
  # variable
  v <- which(values$variables==variable)
  if(length(v)==0) {
    v <- 1
    variable <- values$variable[v]
  }
  
  message('evaluating dvh for ', voi, ', variable: ', variable, '...')
  
  # voi
  index.voi <- get.voi.logical(vois, voi=voi)
  
  # volume
  dx <- mean(diff(values$x))
  dy <- mean(diff(values$y))
  dz <- mean(diff(values$z))
  
  # values
  if(values$Nv>1) {values$values <- values$values[v,,,]}
  voi.values <- values$values[index.voi]
  Nvoxel <- length(voi.values)
  volume.tot <- Nvoxel*dx*dy*dz
  
  # dvh
  voi.values <- sort(voi.values, na.last=FALSE) # mette i valori NA davanti
  volume=sort((1:Nvoxel)/Nvoxel, decreasing=TRUE)
  # elimina i valori NA
  if(!with.na) {
    volume <- volume[!is.na(voi.values)]
    voi.values <- voi.values[!is.na(voi.values)]
  }
  if(with.0) {
    #volume <- c(1, volume)
    if(is.numeric(with.0)) {message('using with.0: ', with.0); min.voi.values <- with.0} else {min.voi.values <- 0}
    volume <- c(max(volume), volume)
    voi.values <- c(min.voi.values, voi.values)
  }
  
  # sovrascrive il nome del voi
  if(is.null(name)) {name <- voi}
    
  return(list(value=voi.values, volume=volume, volume.tot=volume.tot, variable=variable, voi=name, Nvoxel=Nvoxel))
}


#' calcola DVH per tutti i VOI
#' 
#' OBSOLETA ... da integrare con dvh.evaluate
#' 
#' @family DVH
#' @export
dvh.evaluate.all <- function(values=NULL, vois=NULL, variable=NULL, with.0=FALSE, with.na=FALSE)
{
  dvh.all <- list()
  for(v in 1:length(vois$vois)) {
    dvh.all[[v]] <- dvh.evaluate(values=values,
                                 vois=vois,
                                 variable=variable,
                                 voi=vois$vois[v],
                                 with.0=with.0,
                                 with.na=with.na)
  }
  return(dvh.all)
}


#' Evaluate DVH confidence band
#' 
#' Evaluate the DVH confidence band of a distribution (a list) of DVHs, for a specified alpha (percentile) value.
#' 
#' @param dvh.list list of DVHs
#' @param alpha the alpha value
#' @return a list of DVHs consisting of:
#' \item{dvh.mean}{mean}
#' \item{dvh.median}{median}
#' \item{dvh.alpha.lo}{lower band}
#' \item{dvh.alpha.up}{upper band}
#' \item{dvh.min}{min values}
#' \item{dvh.max}{max values}
#' 
#' @family DVH
#' @export
dvh.evaluate.bands <- function(dvh.list, alpha=0.37) {
  Nd <- length(dvh.list)
  v.list <- list()
  # crea lista "semplice di vettori"
  for(i in 1:Nd) {
    v.list[[i]] <- dvh.list[[i]]$value # calcola sui valori, non sui volumi
  }
  v.band <- evaluate.bands(v.list, alpha=alpha)
  message('band -> alpha=', alpha)

  # crea nuovi dvh prendendo come template i vecchi
  dvh.mean <-
  dvh.min <- dvh.max <-
  dvh.median <- dvh.alpha.lo <-
  dvh.alpha.up <- dvh.list[[1]]
  dvh.mean$value <- v.band$mean
  dvh.max$value <- v.band$max
  dvh.min$value <- v.band$min
  dvh.max$value <- v.band$max
  dvh.median$value <- v.band$median
  dvh.alpha.lo$value <- v.band$alpha.lo
  dvh.alpha.up$value <- v.band$alpha.up
  
  return(list(dvh.mean=dvh.mean, dvh.min=dvh.min, dvh.max=dvh.max,
              dvh.median=dvh.median, dvh.alpha.lo=dvh.alpha.lo,
              dvh.alpha.up=dvh.alpha.up))
}


# UTILITIES --------------------------------------------------------------------


#' Evaluate percentile band
#' 
#' Evaluate the percentile (confidence) band for a list of vectors.
#' 
#' @param v.list the list of vectors. Note: the vectors should have the same number of elements.
#' @param alpha the alpha value.
#' @return a list of vectors, consisting of:
#' \item{mean}{the mean values}
#' \item{median}{the median values}
#' \item{alpha.lo}{lower confidence band}
#' \item{alpha.up}{upper confidence band}
#' \item{min}{min values}
#' \item{max}{max values}
#' 
#' @family Utilities
#' @export
evaluate.bands <- function(v.list, alpha=0)
{
  Nv <- length(v.list)
  if(Nv<2) {
    message('number of vectors < 2: no evaluation of bands.')
    return(v.list)
  }
  
  # crea array equivalente
  for(i in 1:Nv) {
    a.tmp <- v.list[[i]]
    if(i==1) {a <- a.tmp} else {a <- rbind(a, a.tmp)}
  }
  
  v.q <- apply(a,2, function(x) quantile(x, c(alpha, 0.5, 1-alpha)))
  v.mean <- apply(a, 2, FUN='mean')
  v.max <- apply(a, 2, FUN='max')
  v.min <- apply(a, 2, FUN='min')
  
  if(alpha==0) {
    message('using alpha=0 in bands...')
    v.q[1,] <- v.min
    v.q[3,] <- v.max
  }  
  
  return(list(mean=v.mean, median=v.q[2,], alpha.lo=v.q[1,], alpha.up=v.q[3,], max=v.max, min=v.min))
}

