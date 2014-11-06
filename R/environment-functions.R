# ENVIRONMENT ------------------------------------------------------------------

# setta le variabili visibili a tutte le funzioni (da http://www.r-bloggers.com/package-wide-variablescache-in-r-packages/)
# sembra che il codice fuori dalle funzioni sia automaticamente eseguito quando si carica il pacchetto
# nota potrei mettere questa inizializzazione in un file setenv.R da mettere nella cartella di lavoro.
# Appena caricate le librerie, potrebbe eseguire automaticamente questo script per settare le variabili.
# se per caso non trova questo file da solo un warning e usa un settaggio predefinito.
# potrebbe esistere una funzione per generare questo file localemente, con i settaggi correnti, in maniera da poterlo
# editare al bisogno.
dektoolsEnv <- new.env()

#' Set the system environment variables
#' 
#' Set the esystem environment variables for Plankit, Gate and Survival.
#' @family Environment
#' @export
setenv.rplanit <- function()
{
  # HOME
  #my.home <- paste(Sys.getenv('HOME'), 'R', sep='/')
  my.home <- get.install.path()
  
  # Gate
  Sys.setenv(DYLD_LIBRARY_PATH=paste(my.home, 'Gate6.2-install/root_v5.34/lib', sep='/'),
             G4ABLADATA=paste(my.home, 'Gate6.2-install/geant4.9.5.p01-install/share/Geant4-9.5.1/data/G4ABLA3.0', sep='/'),
             G4LEDATA=paste(my.home, 'Gate6.2-install/geant4.9.5.p01-install/share/Geant4-9.5.1/data/G4EMLOW6.23', sep='/'),
             G4LEVELGAMMADATA=paste(my.home, 'Gate6.2-install/geant4.9.5.p01-install/share/Geant4-9.5.1/data/PhotonEvaporation2.2', sep='/'),
             G4NEUTRONHPDATA=paste(my.home, 'Gate6.2-install/geant4.9.5.p01-install/share/Geant4-9.5.1/data/G4NDL4.0', sep='/'),
             G4NEUTRONXSDATA=paste(my.home, 'Gate6.2-install/geant4.9.5.p01-install/share/Geant4-9.5.1/data/G4NEUTRONXS1.1', sep='/'),
             G4PIIDATA=paste(my.home, 'Gate6.2-install/geant4.9.5.p01-install/share/Geant4-9.5.1/data/G4PII1.3', sep='/'),
             G4RADIOACTIVEDATA=paste(my.home, 'Gate6.2-install/geant4.9.5.p01-install/share/Geant4-9.5.1/data/RadioactiveDecay3.4', sep='/'),
             G4REALSURFACEDATA=paste(my.home, 'Gate6.2-install/geant4.9.5.p01-install/share/Geant4-9.5.1/data/RealSurface1.0', sep='/'),
             LIBPATH=paste(my.home, 'Gate6.2-install/root_v5.34/lib', sep='/'),
             PYTHONPATH=paste(my.home, 'Gate6.2-install/root_v5.34/lib', sep='/'),
             ROOTSYS=paste(my.home, 'Gate6.2-install/root_v5.34', sep='/'),
             SHLIB_PATH=paste(my.home, 'Gate6.2-install/root_v5.34/lib', sep='/')
  )
  
  ld_library_path <- paste(Sys.getenv('LD_LIBRARY_PATH'),
                           ':', my.home, '/Gate6.2-install/geant4.9.5.p01-install/lib',
                           ':', my.home, '/Gate6.2-install/root_v5.34/lib',
                           ':', my.home, '/Gate6.2-install/2.1.1.0/CLHEP/lib', sep='')
  manpath <- paste(Sys.getenv('MANPATH'),
                   ':', my.home, '/Gate6.2-install/root_v5.34/man', sep='')
  path <- paste(Sys.getenv('PATH'),
                ':', my.home, '/Gate6.2-install/geant4.9.5.p01-install/bin',
                ':', my.home, '/Gate6.2-install/root_v5.34/bin',
                ':', my.home, '/Gate6.2-install/2.1.1.0/CLHEP/bin',
                ':', my.home, '/Gate6.2-install/gate.6.2-install/bin',
                ':', my.home, '/Gate6.2-install/sanitize-hdr', sep='')
  Sys.setenv(LD_LIBRARY_PATH=ld_library_path, MANPATH=manpath, PATH=path)
  
  # Pure-dek
  Sys.setenv(LOG4CXX_CONFIGURATION=paste0(my.home, '/Pure-dek-install/log4cxx.properties'),
             PATH_TO_WEPL=paste0(my.home, '/Pure-dek-install/Beams/'),
             PATH_TO_LUT=paste0(my.home, '/Pure-dek-install/LUT/'))
  path <- paste0(Sys.getenv('PATH'), ':', my.home, '/Pure-dek-install')
  Sys.setenv(PATH=path)
  libs <- paste0(Sys.getenv('LD_LIBRARY_PATH'), ':', my.home, '/Pure-dek-install/ext-libs')
  Sys.setenv(LD_LIBRARY_PATH=libs) # Linux
  dyld.libs <- paste0(Sys.getenv('DYLD_LIBRARY_PATH'), ':', my.home, '/Pure-dek-install/ext-libs')
  Sys.setenv(DYLD_LIBRARY_PATH=dyld.libs) # mac...
  #dyld.libs <- paste0(Sys.getenv('DYLD_FALLBACK_LIBRARY_PATH'), ':', my.home, '/Pure-dek-install/ext-libs')
  #Sys.setenv(DYLD_FALLBACK_LIBRARY_PATH=dyld.libs) # mac (alternativa)...
  
  # Survival
  Sys.setenv(DATA=paste0(my.home, '/Survival-install/data/'))
  path <- paste0(Sys.getenv('PATH'), ':', my.home, '/Survival-install')
  ld_library_path <- paste0(Sys.getenv('LD_LIBRARY_PATH'), ':', my.home, '/Survival-install/lib')
  Sys.setenv(PATH=path, LD_LIBRARY_PATH=ld_library_path)
  
  # NAMESPACE
  assign('dek.setenv', paste0(my.home, '/DEK-install/setenv.sh'), envir=dektoolsEnv)
  #assign('lem.setenv', paste0(my.home, '/Survival/setenv.sh'), envir=dektoolsEnv) 
  assign('lem.setenv', paste0(my.home, '/Survival-install/setenv.sh'), envir=dektoolsEnv)
  assign('gate.setenv', paste0(my.home, '/Gate6.2-install/setenv.sh'), envir=dektoolsEnv)
  assign('gate.template', paste0(my.home, '/Gate6.2-install/tps.template'), envir=dektoolsEnv)
}

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


#'Get installing folder
#'
#'Get the path to the folder where external precompiled packages (such as plankit, survivals and gate) are/will be installed.
#'@export
#'@family Environmen
get.install.path <- function()
{
  # usa la variabile d'ambiente $R_LIBS_USER
  return(paste0(path.expand(Sys.getenv('R_LIBS_USER')), '/'))
}

#'Check the configuration of PlanKIT
#'
#'If the check fails it generates an error.
#'@param use.warning if TRUE it uses warning instead of an error.
#'@export
#'@family Environment
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
#'@family Environment
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

#' Set LUT folder
#' 
#' Set the LUT (look up table) folder to be used in pure-dek computations.
#' 
#' @param lut.folder The lut.folder. If it is not defined the lut folder is set to the default one (the pure-dek local installation.)
#' @export
#' @family Environment
set.lut.folder <- function(lut.folder=NULL)
{
  if(is.null(lut.folder)) {
    my.home <- get.install.path()
    lut.folder <- paste0(my.home, '/Pure-dek-install/LUT/')
  }
  
  Sys.setenv(PATH_TO_LUT=path.expand(lut.folder))
}

# INSTALLING -------------------------------------------------------------------

#' Install Gate
#' @export
#' @family Install
install.gate <- function()
{
  # installing path
  my.home <- get.install.path()
  dir.create(my.home, recursive=TRUE, showWarnings=FALSE)
  
  # check del sistema operativo
  my.system <- Sys.info()['sysname']
  if(my.system=='Linux') {
    message('installing compiled plakit for Linux...')
  } else {
    stop(paste0('error: there is no precompiled available for your system (', my.system, ')'))
  }
  
  message('downloading gate...')
  download.file(url='http://totlxl.to.infn.it/tools/Gate6.2-install.tar.bz2', destfile='Gate6.2-install.tar.bz2')
  
  message('uncompressing gate...')
  system('tar jxf Gate6.2-install.tar.bz2; rm Gate6.2-install.tar.bz2', ignore.stdout=TRUE, ignore.stderr=TRUE)
  
  install.dir <- paste(my.home, 'Gate6.2-install', sep='/')
  message('moving gate in ', install.dir, ' ...')
  unlink(install.dir, recursive=TRUE)
  file.rename(from='Gate6.2-install', to=install.dir)
}


#' Install Plankit (obsolete)
#' @export
#' @family Install
install.plankit <- function()
{stop('install.plankit() is obsolete... Use install.puredek() instead ;)')}

#' Install Pure-dek
#' @export
#' @family Install
install.puredek <- function()
{
  # installing path
  my.home <- get.install.path()
  dir.create(my.home, recursive=TRUE, showWarnings=FALSE)
  
  # check del sistema operativo
  my.system <- Sys.info()['sysname']
  if(my.system=='Linux') {
    message('installing compiled Pure-dek for Linux...')
    puredek_package <- 'Pure-dek-install.tar.bz2'
  } else if(my.system=='Darwin') {
    message('installing compiled Pure-dek for Mac (Darwin)...')
    puredek_package <- 'Pure-dek-install.mac.tar.bz2'
  } else {
    stop(paste0('error: there is no precompiled available for your system (', my.system, ')'))
  }
  
  message('downloading pure-dek...')
  download.file(url=paste0('http://totlxl.to.infn.it/tools/', puredek_package), destfile='Pure-dek-install.tar.bz2')
  
  message('uncompressing pure-dek...')
  system('tar jxf Pure-dek-install.tar.bz2; rm Pure-dek-install.tar.bz2', ignore.stdout=TRUE, ignore.stderr=TRUE)
  
  install.dir <- paste(my.home, 'Pure-dek-install', sep='/')
  message('moving pure-dek in ', install.dir, ' ...')
  unlink(install.dir, recursive=TRUE)
  file.rename(from='Pure-dek-install', to=install.dir)
}

#' Install Survival
#' @export
#' @family Install
install.survival <- function()
{
  # installing path
  my.home <- get.install.path()
  dir.create(my.home, recursive=TRUE, showWarnings=FALSE)
  
  # check del sistema operativo
  my.system <- Sys.info()['sysname']
  if(my.system=='Linux') {
    message('installing compiled survival for Linux...')
    survival_package <- 'Survival-install.tar'
  } else if(my.system=='Darwin') {
    message('installing compiled survival for Mac (Darwin)...')
    survival_package <- 'Survival-install.mac.tar'
  } else {
    stop(paste0('error: there is no precompiled available for your system (', my.system, ')'))
  }
  
  message('downloading survival...')
  download.file(url=paste0('http://totlxl.to.infn.it/tools/', survival_package, '.bz2'), destfile='Survival-install.tar.bz2')
  
  message('uncompressing survival...')
  system(paste0('tar jxf Survival-install.tar.bz2; rm Survival-install.tar.bz2'), ignore.stdout=TRUE, ignore.stderr=TRUE)
  
  install.dir <- paste(my.home, 'Survival-install', sep='/')
  message('moving survival in ', install.dir, ' ...')
  unlink(install.dir, recursive=TRUE)
  file.rename(from='Survival-install', to=install.dir)
}

# RUNNING ----------------------------------------------------------------------
# qui vengono eseguite tutte procedure al momento del caricamento del pacchetto.

.onLoad <- function(libname=find.package("Rplanit"), pkgname = "Rplanit")
{
  setenv.rplanit()
}
