# TESTING FUNCTIONS ------------------------------------------------------------

#' carica dati per i test
#' @family Testing Functions
load.test.data <- function()
{
  # esempio
  plan <- read.plan('data-test/waterbox-plan')
  values <- get.values(plan)
  vois <- get.vois(plan)
  contours <- get.contours(plan)
  ct <- get.ct(plan)
  
  # ct + contorni "reali"
  ct <- read.3d.array('../../analisi/CTs-DICOM/HN/plankit/ds.3d')
  contours <- read.contours(file.contours='../../analisi/CTs-DICOM/HN/plankit/ds_R3PIDEd0.2v0+MKM.contours', file.CT='../../analisi/CTs-DICOM/HN/plankit/ds.3d')
  
  return(list(plan, values, vois, contours))
}

# test per generare la CT
ct.generate.test <- function() {
  CT.df <- data.frame(x_min=-25, x_max=+25, y_min=-25, y_max=+25, z_min=-25, z_max=+25, HU=-100, voi='BODY')
  CT.df <- rbind(CT.df,
                 data.frame(x_min=-5, x_max=+5, y_min=-5, y_max=+5, z_min=-5, z_max=+5, HU=0, voi='PTV'))
  CT <- generate.ct(CT.df=CT.df, deltaX=1, deltaY=1, deltaZ=1)
  display.slice.ct(CT)
}

# carica direttamente le funzioni
load.dektools <- function()
{
  # librerie esterne
  library(knitr)
  library(ggplot2)
  library(fields)
  library(splancs)
  library(markdown)
  library(rgl)
  library(misc3d)
  
  # carica tutte le funzioni definite nel folder...
  
  # variabili di environment
  my.folder  <- '/home/andrea/Programmi/dektools/R/'
  dek.setenv <- '/home/andrea/dek/trunk/DEK/bin/setenv.sh'
  lem.setenv <- '/opt/dek-tools/lem-setenv.sh'
  gate.setenv <- '/home/andrea/Gate6.2/setenv.sh'
  gate.template <- '/home/andrea/Gate6.2/tps.template'
  
  my.files <- dir(my.folder, pattern=glob2rx('*.R'))
  my.files <- my.files[my.files!='load-functions.R']
  
  for(my.file in my.files) {
    source(paste(my.folder, my.file, sep=''))
  }
}
