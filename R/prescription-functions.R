# PRESCRIPTION -----------------------------------------------------------------


#' Crea oggetto prescrizione (dataframe)
#' 
#' nota: si tratta di un dataframe: diversi fields sono semplicemente aggiunti
#' appendendo delle righe al dataframe...
#' 
#' @family Prescription
#' @export
#' 
create.constraint <- function(N=1,
                              VOI='PTV',
                              VOIIndex=NA,
                              type='EXACT',
                              variable='Dose[Gy]',
                              value=2.00,
                              volumeFraction=1.0,
                              weight=1.0)
{
  constraint <- data.frame(constraint=1:N,
                           VOI=VOI,
                           VOIIndex=VOIIndex,
                           type=type,
                           variable=variable,
                           value=value,
                           volumeFraction=volumeFraction,
                           weight=weight,
                           stringsAsFactors=FALSE)
  return(constraint)
}


#' recupera prescrizione
#' 
#' @family Prescription
#' @export
#' 
get.prescription <- function(plan)
{
  return(plan$prescription)
}


#' Check the prescription
#' 
#' It verifies if the contraints defined in the prescription are verified. It uses the actual values distributions to perform the check.
#' 
#' @param plan the plan object
#' @param values the values object. If it is \code{NULL}, \code{values} is obtained directly from the \code{plan}
#' @param vois the vois object. If it is \code{NULL}, \code{vois} is obtained directly from the \code{plan}
#' @param prescription the prescription dataframe. If it is \code{NULL}, \code{prescription} is obtained directly from the \code{plan}
#' 
#' @family Prescription
#' @export

check.prescription <- function(plan=plan, values=NULL, vois=NULL, prescription=NULL)
{
  # controllo degli oggetti passati
  if(is.null(prescription)) {prescription <- get.prescription(plan)}
  if(is.null(vois)) {vois <- get.vois(plan)}
  if(is.null(values)) {values <- get.values(plan)}
  
  check <- rep(NA, nrow(prescription))
  for(i in 1:nrow(prescription)){
    message('checking prescription for voi: ', prescription$VOI[i])
    dvh <- dvh.evaluate(values=values, vois=vois, voi=prescription$VOI[i], variable=prescription$variable[i])
    if(prescription$type[i]=='EXACT') {
      check[i] <- sum(dvh$value==prescription$value[i], na.rm=TRUE)/sum(!is.na(dvh$value)) == prescription$volumeFraction[i]
    }
    if(prescription$type[i]=='UPPER') {
      check[i] <- max(dvh$volume[dvh$value >= prescription$value[i]]) <= prescription$volumeFraction[i]
    }
    if(prescription$type[i]=='LOWER') {
      check[i] <- min(dvh$volume[dvh$value <= prescription$value[i]]) >= prescription$volumeFraction[i]
    }
  }
  return(cbind(prescription, check))
}
