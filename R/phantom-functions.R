#' Generate a CT phantom
#' 
#' It uses the definition of box-shaped volumes stored in a data.frame with columns:
#' x_min, x_max, y_min, y_max, z_min, z_max, HU, voi.
#' 
#' @param CT.df the dataframe in which are defined the volumes.
#' @param deltaX, deltaY, deltaZ are the dimensions of the voxels.
#' 
#' @return a CT (values) object
#' 
#' @export

generate.ct <- function(CT.df, deltaX, deltaY, deltaZ)
{
  numberOfVolumes <- nrow(CT.df)
  
  # calcola numero di voxels
  message('evaluating number of voxels...')
  x.min <- min(CT.df$x_min)
  x.max <- max(CT.df$x_max)
  y.min <- min(CT.df$y_min)
  y.max <- max(CT.df$y_max)
  z.min <- min(CT.df$z_min)
  z.max <- max(CT.df$z_max)
  numberOfX <- round((x.max-x.min)/deltaX)
  numberOfY <- round((y.max-y.min)/deltaY)
  numberOfZ <- round((z.max-z.min)/deltaZ)
  numberOfVoxels <- numberOfZ * numberOfX * numberOfY
  message('creating array ', numberOfX, 'x', numberOfY, 'x', numberOfZ, ' = ', numberOfVoxels)
  
  # calcola coordinate voxel "zero"
  x0 <- min(CT.df$x_min) + deltaX/2
  y0 <- min(CT.df$y_min) + deltaY/2
  z0 <- min(CT.df$z_min) + deltaZ/2
  
  # crea CT d'acqua
  x <- (0:(numberOfX-1))*deltaX + x0
  y <- (0:(numberOfY-1))*deltaY + y0
  z <- (0:(numberOfZ-1))*deltaZ + z0
  CT <- array(rep(0, numberOfVoxels), dim=c(numberOfX, numberOfY, numberOfZ))
  
  # riempie i sottovolumi
  for (i in 1:numberOfVolumes) {
    i_min <- ceiling((CT.df$x_min[i] - x0)/deltaX) + 1
    i_max <- floor((CT.df$x_max[i] - x0)/deltaX) + 1
    j_min <- ceiling((CT.df$y_min[i] - y0)/deltaY) + 1
    j_max <- floor((CT.df$y_max[i] - y0)/deltaY) + 1
    k_min <- ceiling((CT.df$z_min[i] - z0)/deltaZ) + 1
    k_max <- floor((CT.df$z_max[i] - z0)/deltaZ) + 1
    
    #print(c(i_min, i_max, j_min, j_max, k_min, k_max))
    
    CT[i_min:i_max, j_min:j_max, k_min:k_max] <-CT.df$HU[i]
  }
  
  # crea struttura CT
  CT.values <- list(values=CT, x=x, y=y, z=z, Nx=numberOfX, Ny=numberOfY, Nz=numberOfZ, Nv=1, variables='HounsfieldNumber')
  
  # ritorna la struttura
  return(CT.values)
}


#' Generate the contours
#' 
#' It generates a contours object from the definition of the VOIs.
#' The VOIs are box-shaped an are defined in a data.frame with columns:
#' x_min, x_max, y_min, y_max, z_min, z_max, voi. Each row correspond to a VOI.
#' It needs, to identify the slices (axial slices), a reference CT. Alternatively the z vector can be used directly. Optionally a tissue name can be supplied.
#' 
#' @param contours.df the dataframe in which the definition of the volumes is stored.
#' @param CT a CT object.
#' @param z a vector that specifies the z of the axial planes.
#' @param tissue the tissue name.
#' 
#' @return a contours object
#' 
#' @export

generate.contours <- function(contours.df, CT=NULL, z=NULL, tissue='R3PIDEd0.2v0+MKM')
{
  numberOfvois <- nrow(contours.df)
  
  if(!is.null(CT)) {
    numberOfZ <- CT$Nz
    z <- CT$z
  } else {
    numberOfZ <- length(z)
  }
  
  vois <- paste(contours.df$voi, '_', tissue, sep='')
  
  for(i in 1:numberOfvois) {
    id.polygon <- 0
    for(iz in 1:numberOfZ) {
      if(z[iz]>=contours.df$z_min[i] & z[iz]<=contours.df$z_max[i]) {
        xc <- c(contours.df$x_min[i], contours.df$x_max[i], contours.df$x_max[i], contours.df$x_min[i])
        yc <- c(contours.df$y_min[i], contours.df$y_min[i], contours.df$y_max[i], contours.df$y_max[i])
        contours.tmp <- data.frame(id=(i-1), polygon=id.polygon, slice=iz-1, x=xc, y=yc, z=z[iz], contour=contours.df$voi[i], tissue=tissue, type=contours.df$voi[i])
        id.polygon <- id.polygon+1
        if(i==1 & id.polygon==1) {contours <- contours.tmp} else {contours <- rbind(contours, contours.tmp)}
      }
    }
  }
  
  return(contours)
}


#' Generate a complete phantom (CT + contours)
#' 
#' It generates a CT and assocuated contours, via the definition of box-shaped volumes.
#' 
#' It uses the definition of volumes and VOIs stored in two dataframes:
#' 
#' CT.df:
#' x_min, x_max, y_min, y_max, z_min, z_max, HU, voi.
#' 
#' contours.df:
#' x_min, x_max, y_min, y_max, z_min, z_max, voi.
#' 
#' It is possible to use directly the same CT.df for the definition of the volumes of the contours (contours.df).
#' 
#'CT and contours are directly written in two files with names specfied by \code{name}.
#' 
#' @param the dataframe in which the definition of the volumes for the CT is stored
#' @param the dataframe in which the definition of the volumes for the contours is stored
#' @param deltaX, deltaY, deltaZ are the dimensions of the voxels
#' @param the name of the phantom. Files are automaticallywritten for the CT (<name>.3d) and the contours (<name>.contours)
#' @param tissue the tissue name
#' 
#' @export
#' 
generate.ct.contours <- function(CT.df, contours.df, deltaX, deltaY, deltaZ, name='waterbox', tissue='R3PIDEd0.2v0+MKM')
{
  # crea struttura CT e la salva
  CT.values <- generate.ct(CT.df, deltaX, deltaY, deltaZ)
  write.3d.array(CT.values, paste(name, '.3d', sep=''))

  # Calcola i contorni e li salva
  contours <- generate.contours(CT=CT.values, contours.df=contours.df, tissue=tissue)
  write.contours(contours=contours, name=name)

}
