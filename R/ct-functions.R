# CT ---------------------------------------------------------------------------

#' Get ct object from plan object
#' @param plan The plan object (it can be a list of plans)
#' @return A CT object (or a list of CT objects)
#' @family CT
#' @export
get.ct <- function(plan) UseMethod("get.ct")
{
  #   if(is.null(plan[['outputValuesFile']])) {
  #     cat('The plan "', plan[['name']], '" has no values file.\n', sep='')
  #     return(NULL)
  #   }
  
  #return(read.3d(plan$ctFile))
}

#' @family CT
#' @export
get.ct.list <- function(plans)
{
  CTs <- list()
  for(i in 1:length(plans)) {
    CTs[[i]] <- get.ct(plans[[i]])
  }
  return(CTs)
}

#' @family CT
#' @export
get.ct.plankit.plan <- function(plan)
{
  if(!is.null(plan[['ct']])) {
    return(plan[['ct']])
  } else {
    return(read.3d(plan$ctFile))
  }
}

#' @family CT
#' @export
get.ct.gate.plan <- function(plan)
{
  if(!is.null((plan[['ct']]))) {return(plan[['ct']])}
  else if(!is.null(plan$ctFile)) {return(read.3d.hdr(file.name=plan$ctFile))}
  else if(!is.null(plan$ctWater)) {
    CT.df <- data.frame(x_min=-plan$ctWater$Dx/2, x_max=plan$ctWater$Dx/2,
                        y_min=-plan$ctWater$Dy/2, y_max=plan$ctWater$Dy/2,
                        z_min=-plan$ctWater$Dz/2, z_max=plan$ctWater$Dz/2,
                        HU=1, voi='BODY')
    ct <- generate.ct(CT.df=CT.df, deltaX=plan$voxelSize[1], deltaY=plan$voxelSize[2], deltaZ=plan$voxelSize[3])
    return(ct)
  }
  stop('No method for class plan.gate')
}


#' rimuove valori fuori range della ct
#' 
#' @family CT
#' @export
sanitize.ct <- function(ct, HU.min=-1000, HU.max=3000, threshold.air=NULL)
{
  ct$values[ct$values<HU.min] <- HU.min
  ct$values[ct$values>HU.max] <- HU.max
  if(!is.null(threshold.air)) {ct$values[ct$values<=threshold.air] <- HU.min}
  return(ct)
}
