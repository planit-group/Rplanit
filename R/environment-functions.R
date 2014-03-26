# ENVIRONMENT ------------------------------------------------------------------

# setta le variabili visibili a tutte le funzioni (da http://www.r-bloggers.com/package-wide-variablescache-in-r-packages/)
# sembra che il codice fuori dalle funzioni sia automaticamente eseguito quando si carica il pacchetto
# nota potrei mettere questa inizializzazione in un file setenv.R da mettere nella cartella di lavoro.
# Appena caricate le librerie, potrebbe eseguire automaticamente questo script per settare le variabili.
# se per caso non trova questo file da solo un warning e usa un settaggio predefinito.
# potrebbe esistere una funzione per generare questo file localemente, con i settaggi correnti, in maniera da poterlo
# editare al bisogno.
dektoolsEnv <- new.env()
assign('dek.setenv', '/home/andrea/dek/trunk/DEK/bin/setenv.sh', envir=dektoolsEnv)
assign('lem.setenv', '/opt/dek-tools/lem-setenv.sh', envir=dektoolsEnv)
assign('gate.setenv', '/home/andrea/Gate6.2/setenv.sh', envir=dektoolsEnv)
assign('gate.template', '/home/andrea/Gate6.2/tps.template', envir=dektoolsEnv)

#' Set Environment Variables
#' 
#' @param dek.setenv path to configuration script plankit TPS.
#' @param lem.setenv path to the cofiguration script for the radiobiological simulation code.
#' @param gate.setenv path to the configuration script for Gate/Geant4 code.
#' @param gate.template path to the folder containing the template macro files for the Gate4/Geant4 simulations.
#' @export
#' @family Environment
set.dektools.env <- function(dek.setenv=NULL, lem.setenv=NULL, gate.setenv=NULL, gate.template=NULL)
{
  if(!is.null(dek.setenv)) {assign('dek.setenv', dek.setenv, envir=dektoolsEnv)}
  if(!is.null(lem.setenv)) {assign('lem.setenv', lem.setenv, envir=dektoolsEnv)}
  if(!is.null(gate.setenv)) {assign('gate.setenv', gate.setenv, envir=dektoolsEnv)}
  if(!is.null(gate.template)) {assign('gate.template', gate.template, envir=dektoolsEnv)}
  
}

#' Get Environment Variables
#' 
#' @export
#' @family Environment
get.dektools.env <- function()
{
  return(list(dek.setenv=get('dek.setenv', envir=dektoolsEnv),
              lem.setenv=get('lem.setenv', envir=dektoolsEnv),
              gate.setenv=get('gate.setenv', envir=dektoolsEnv),
              gate.template=get('gate.template', envir=dektoolsEnv)))
}


#'Check the configuration of PlanKIT
#'
#'If the check fails it generates an error.
#'@param use.warning if TRUE it uses warning instead of an error.
#'@export
#'@faily Environment
check.plankit <- function(use.warning=FALSE)
{
  dek.setenv <- get('dek.setenv', envir=dektoolsEnv)
  if(file.exists(dek.setenv)) {
    message(paste('Variable dek.setenv correctly set (\"', dek.setenv, '\")\n', sep=''))
  } else {
    if(use.warning) {
      warning('Error in variable dek.setenv (\"', dek.setenv, '\" not found)')
    } else {
      stop('Error in variable dek.setenv (\"', dek.setenv, '\" not found)')
    }
  }
}

#'Check the configuration of Gate
#'
#'If the check fails it generates an error.
#'#'@param use.warning if TRUE it uses warning instead of an error.
#'@export
#'@faily Environment
check.gate <- function(use.warning=FALSE)
{
  gate.setenv <- get('gate.setenv', envir=dektoolsEnv)
  gate.template <- get('gate.template', envir=dektoolsEnv)
  if(file.exists(gate.setenv)) {
    message(paste('Variable gate.setenv correctly set (\"', gate.setenv, '\")\n', sep=''))
  } else {
    if(use.warning) {
      warning('Error in variable gate.setenv (\"', gate.setenv, '\" not found)')
    } else {
      stop('Error in variable gate.setenv (\"', gate.setenv, '\" not found)')
    }
  }
  if(file.exists(gate.template)) {
    message(paste('Variable gate.template correctly set (\"', gate.template, '\")\n', sep=''))
  } else {
    if(use.warning) {
      warning('Error in variable gate.template (\"', gate.template, '\" not found)')
    } else {
      stop('Error in variable gate.template (\"', gate.template, '\" not found)')
    }
  }
}