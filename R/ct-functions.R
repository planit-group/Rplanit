# CT ---------------------------------------------------------------------------

#' recupera ct
#' 
#' @family CT
#' @export
get.ct <- function(plan)
{
  #   if(is.null(plan[['outputValuesFile']])) {
  #     cat('The plan "', plan[['name']], '" has no values file.\n', sep='')
  #     return(NULL)
  #   } 
  return(read.3d(plan$ctFile))
}


#' rimuove valori fuori range della ct
#' 
#' @family CT
#' @export
sanitize.ct <- function(ct, HU.min=-1000, HU.max=3000)
{
  ct$values[ct$values<HU.min] <- HU.min
  ct$values[ct$values>HU.max] <- HU.max
  return(ct)
}
