# TESTING FUNCTIONS ------------------------------------------------------------

#' Evaluate a test planning
#' 
#' Generate a simple water phantom CT with a central PTV and four OARs around.
#' Set a plan with protons (default) with 4 coplanar fields and perform an optimization to deliver exatly 1 Gy(RBE) (default) in the PTV.
#' Generate a graphical report.
#' 
#' @param plan.name The name of the plan.
#' @param particle The primary particle ('proton', 'lithium', 'helium', 'carbon' or 'oxygen')
#' @param dose The equivalent dose to be delivered in the PTV (Gy(RBE))
#' @param outmessages Display detailed output messages of the evaluation. Warning: the throughput is high, use this option only when running from a simple text terminal.
#' @return The plan object.
#' @export
#' @family Testing Functions
demo.puredek <- function(plan.name='test-plan', particle='proton', dose=1, outmessages=FALSE) {
  
  # crea CT
  CT.df <- data.frame(x_min=-85, x_max=85, y_min=-85, y_max=85, z_min=-85, z_max=85, HU=0)
  ct <- generate.ct(CT.df=CT.df, deltaX=1, deltaY=1, deltaZ=1)
  
  # crea contorni
  contours.BODY <- data.frame(x_min=-85, x_max=85, y_min=-85, y_max=85, z_min=-85, z_max=85, voi='BODY')
  contours.PTV <- data.frame(x_min=-15, x_max=15, y_min=-15, y_max=15, z_min=-15, z_max=15, voi='PTV')
  
  dd <- 35
  contours.OAR1 <- data.frame(x_min=-15-dd, x_max=15-dd, y_min=-15, y_max=15, z_min=-15, z_max=15, voi='OAR1')
  contours.OAR2 <- data.frame(x_min=-15+dd, x_max=15+dd, y_min=-15, y_max=15, z_min=-15, z_max=15, voi='OAR2')
  contours.OAR3 <- data.frame(x_min=-15, x_max=15, y_min=-15-dd, y_max=15-dd, z_min=-15, z_max=15, voi='OAR3')
  contours.OAR4 <- data.frame(x_min=-15, x_max=15, y_min=-15+dd, y_max=15+dd, z_min=-15, z_max=15, voi='OAR4')
  contours.df <- rbind(contours.BODY, contours.PTV, contours.OAR1, contours.OAR2, contours.OAR3, contours.OAR4)
  contours.df$tissue='R3PIDEd0.2v0+MKM'
  cc <- generate.contours(contours.df=contours.df, CT=ct)
  
  ## crea piano
  plan <- create.plan(name=plan.name)
  
  # paziente
  plan$ct <- ct
  plan$contours <- cc
  plan$hounsfieldToDensityFile <- paste0(get.install.path(), '/Pure-dek-install/CT/water.HNToDensity.1d')
  plan$hounsfieldToStoichiometryFile <- paste0(get.install.path(), '/Pure-dek-install/CT/water.HNToStoichiometry.1d')
  
  # fields
  beamLines <- c(proton='protonSimple', lithium='lithiumSimple', helium='heliumSimple', carbon='carbonSimple', oxygen='oxygenSimple')
  beamLine <- beamLines[particle]
  fields <- create.field(N = 4, targetVOI = 'PTV', beamLine = beamLine, iecGantryAngle = c(0, 90, 180, 270), spotsExtensionOutsideTarget = 10, targetIsocenter = c(0, 0, 0), interSpotSpacing = c(1.5, 1.5, 1.5))
  plan$fields <- fields
  
  # tipo di calcolo
  plan$computingGridVoxelSizes <- c(3, 3, 3)
  plan$computingGridCoverage <- 0 # BODY
  plan$computingValues <- c('Dose[Gy]', 'MeanLET[keV/um]', 'DoseAveragedLET[keV/um]', 'RBE', 'BiologicalDose[Gy(RBE)]', 'Survival', 'Alpha[Gy^(-1)]', 'Beta[Gy^(-2)]', 'NumberOfLethalEvents')
  plan$cutOffRadius <- 50
  
  # prescrizione/ottimizzazione
  constraint.PTV <- create.constraint(N = 1, VOI = 'PTV', type = 'EXACT', variable = 'BiologicalDose[Gy(RBE)]', value = dose, volumeFraction = 1, weight = 1)
  plan$prescription <- constraint.PTV
  plan$maxOptimizationIterations <- 50
  
  # calcola
  tt <- system.time(plan.out <- run.dek.inverse(plan, outmessages=outmessages)); 
  
  # genera un report
  generate.report(plan.out, N.slice = 1, html=TRUE)
  browseURL('file://report_temp.html')
  
  message('planning computing time:')
  print(tt)
  
  return(plan.out)
}


# DA SISTEMARE -----------------------------------------------------------------

#' carica dati per i test
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
