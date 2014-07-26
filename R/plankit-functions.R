# PLANKIT VARIABLES ------------------------------------------------------------

#' Available beamlines (PlanKIT)
#' 
#' The available beamlines name identifiers to be used with PlanKIT
#' @return A list of the available beamlines
#' @export
#' @family PlanKITVariables
available.beamlines <- function()
{
  return(
  c('protonSimple',
    'protonSimpleGate',
    'carbonSimple',
    'lithiumSimple',
    'heliumSimple',
    'oxygenSimple')
  )
}
  


#' Available variables (PlanKIT)
#' 
#' The available variables name identifiers that can be evaluated with PlanKIT
#' @return A list of the available variables
#' @export
#' @family PlanKITVariables
available.variables <- function()
{
  
  return(
  c('Dose[Gy]',
    'DoseSquared[Gy^2]',
    'MeanLET[keV/um]',
    'DoseAveragedLET[keV/um]',
    'BiologicalDose[Gy(RBE)]',
    'Survival',
    'RBE',
    'NumberOfLethalEvents',
    'Alpha[Gy^(-1)]',
    'Beta[Gy^(-2)]')
  )
}

#' Available constraints (PlanKIT)
#' 
#' The available constraint types name identifiers tha can be used in the optimization procedure of PlanKIT
#' @return A list of the available variables
#' @export
#' @family PlanKITVariables
available.constraint.types <- function()
{
  return(c('EXACT', 'UPPER', 'LOWER'))
}



# PLAN -------------------------------------------------------------------------

#' Generate Plan template (PlanKIT)
#' 
#' Generate a default (empty) template of the Plan to be used with PlanKIT.
#' @return A Plan (plankit.plan object).
#' 
#' @family Plan
#' @export
#' 
create.plan <- function(plan=NA,
                        name='test',
                        ctFile='ct.3d',
                        ct=NULL,
                        contoursFile='plan.contours',
                        contours=NULL,
                        hounsfieldToDensityFile='fluka-schneider.HNToDensity.1d',
                        hounsfieldToStoichiometryFile='fluka-schneider.HNToStoichiometry.1d',
                        computingGridVoxelSizes=c(2,2,2),
                        computingGridCoverage=0,
                        computingGridCoverageBoundary=0,
                        beamLine='protonCnao',
                        computingValues='Dose[Gy]',
                        cutOffRadius=100,
                        cutEnergyFraction=0.01,
                        fields=NULL,
                        prescription=NULL,
                        maxOptimizationIterations=50,
                        outputBeamsFile=NULL,
                        inputBeamsFile=NULL,
                        beams=NULL)
{
  plan <- list(name=name,
               ctFile=ctFile,
               contoursFile=contoursFile,
               hounsfieldToDensityFile=hounsfieldToDensityFile,
               hounsfieldToStoichiometryFile=hounsfieldToStoichiometryFile,
               computingGridVoxelSizes=computingGridVoxelSizes,
               computingGridCoverage=computingGridCoverage,
               computingGridCoverageBoundary=computingGridCoverageBoundary,
               beamLine=beamLine,
               computingValues=computingValues,
               cutOffRadius=cutOffRadius,
               cutEnergyFraction=cutEnergyFraction,
               fields=fields,
               prescription=prescription,
               maxOptimizationIterations=maxOptimizationIterations,
               outputBeamsFile=outputBeamsFile,
               inputBeamsFile=inputBeamsFile,
               beams=beams # include possibilità di includere direttamente un beams object.
  )
  class(plan) <- 'plankit.plan'
  return(plan)
}


#' Save Plan
#' 
#' Save the Plan object on disk. The name of the folder in which the plan data are saved is derived from the plan name
#' 
#' @family Plan
#' @export
#' 
save.plan <- function(plan) {
  
  plan.file <- paste(plan$name, '/plan.Rdata', sep='')
  dir.create(plan$name, recursive=TRUE, showWarnings=FALSE)
  save(plan, file=plan.file, compress='bzip2')
  cat('plan saved in', plan.file, '\n')
  
}


#' Read Plan
#' 
#' Read a Plan object from disk. The origin folder is derived from the plan name.
#' @param name a string name of the plan.
#' @return the Plan object
#' 
#' @family Plan
#' @export
#' 
read.plan <- function(name)
{
  load(paste(name, '/plan.Rdata', sep=''))
  return(plan)
}




# PLANKIT ----------------------------------------------------------------------


#' Inverse planning (PlanKIT)
#' 
#' Perform an inverse planning with PlanKIT.
#' @param plan the Plan object.
#' @param outmessages if TRUE, it displays the messages from PlanKIT output.
#' @return a Plan Object (including data/information of the evaluated output of the inverse planning, e.g. beams, values, vois, etc.)
#' 
#' @family PlanKIT
#' @export
run.dek.inverse <- function(plan, outmessages=FALSE) {
  
  message('running plankit inverse planning...')
  
  # crea folder
  dir.create(plan$name, recursive=TRUE, showWarnings=FALSE)
  
  # rimuove contenuto del folder (il piano viene completamente sovrascritto)
  unlink(paste(plan$name, '/*', sep=''), recursive = FALSE, force = FALSE)
  
  # CT e CONTORNI
  if(!is.null(plan[['ct']])) {
    message('using ct object in plan...')
    write.3d.array(ct=plan$ct, file.name=paste0(plan$name, '/ct.3d'))
    plan$ctFile <- paste0(plan$name, '/ct.3d')
  }
  if(!is.null(plan[['contours']])) {
    write.contours(contours=plan$contours, name=paste0(plan$name, '/plan'))
    plan$contoursFile <- paste0(plan$name, '/plan.contours')
  }
  
  # input files (path assoluti)
  ctFile <- paste(getwd(), '/', plan[['ctFile']], sep='')
  hounsfieldToDensityFile <- paste(getwd(), '/', plan[['hounsfieldToDensityFile']], sep='')
  hounsfieldToStoichiometryFile <- paste(getwd(), '/', plan[['hounsfieldToStoichiometryFile']], sep='')
  contoursFile <- paste(getwd(), '/', plan[['contoursFile']], sep='')
  plan.config <- paste(plan[['name']], '/plan.config', sep='')
  
  # crea file prescrizione
  plan.prescription <- paste(plan$name, '/plan.prescription', sep='')
  for(i in 1:nrow(plan$prescription)) {
    if(is.na(plan$prescription$VOIIndex[i]))
      plan$prescription$VOIIndex[i] <- get.voiindex(voi=plan$prescription$VOI[i],
                                                    file.contours=plan$contoursFile)
  }
  write.table(plan$prescription[c('VOIIndex', 'type', 'variable', 'value', 'volumeFraction', 'weight')],
              file=plan.prescription,
              col.names=FALSE,
              row.names=FALSE,
              quote=FALSE)
  
  # crea file plan.config
  con <- file(plan.config, "w") # open for writing in text mode
  
  writeLines(paste('ctFile = ', ctFile, '\n'), con=con, sep='')
  writeLines(paste('hounsfieldToDensityFile = ', hounsfieldToDensityFile, '\n'), con=con, sep='')
  writeLines(paste('hounsfieldToStoichiometryFile = ', hounsfieldToStoichiometryFile, '\n'), con=con, sep='')
  writeLines(paste('contoursFile = ', contoursFile, '\n'), con=con, sep='')
  p.computingGridVoxelSizes <- paste(plan[['computingGridVoxelSizes']], collapse=' ')
  writeLines(paste('computingGridVoxelSizes = ', p.computingGridVoxelSizes, '\n', collapse=' '), con=con, sep='')
  writeLines(paste('computingGridCoverage = ', plan[['computingGridCoverage']], '\n'), con=con, sep='')
  writeLines(paste('computingGridCoverageBoundary = ', plan[['computingGridCoverageBoundary']], '\n'), con=con, sep='')
  writeLines(paste('beamLine = ', plan[['beamLine']], '\n'), con=con, sep='')
  writeLines(paste('computingValues = ', plan[['computingValues']], '\n'), con=con, sep='')
  writeLines(paste('cutOffRadius = ', plan[['cutOffRadius']], '\n'), con=con, sep='')
  writeLines(paste('cutEnergyFraction = ', plan[['cutEnergyFraction']], '\n'), con=con, sep='')
  
  # fields
  fields <- plan[['fields']]
  for(i in 1:nrow(fields)) {
    if(is.na(fields$targetVOIIndex[i])) {
      fields$targetVOIIndex[i] <- get.voiindex(fields$targetVOI[i],
                                               file.contours=plan$contoursFile)
    }
    writeLines(paste('targetVOIIndex = ', fields[i,'targetVOIIndex'], '\n'), con=con, sep='')
    writeLines(paste('iecGantryAngle = ', fields[i,'iecGantryAngle'], '\n'), con=con, sep='')
    writeLines(paste('iecPatientSupportAngle = ', fields[i,'iecPatientSupportAngle'], '\n'), con=con, sep='')
    writeLines(paste('interSpotSpacing = ', fields[i,'interSpotSpacing.x'], fields[i,'interSpotSpacing.y'], fields[i,'interSpotSpacing.z'], '\n'), con=con, sep='')
    writeLines(paste('spotsExtensionOutsideTarget = ', fields[i,'spotsExtensionOutsideTarget'], '\n'), con=con, sep='')
    if(!is.na(fields[i,'targetIsocenter.x']) &  !is.na(fields[i,'targetIsocenter.y']) & !is.na(fields[i,'targetIsocenter.z'])) {
      writeLines(paste('isocenter =', fields[i,'targetIsocenter.x'], fields[i,'targetIsocenter.y'], fields[i,'targetIsocenter.z'], '\n'), con=con, sep=' ')
    }
  }
  plan$fields <- fields
  
  # optimization
  writeLines('prescriptionFile = plan.prescription\n', con=con, sep='')
  writeLines(paste('maxOptimizationIterations = ', plan[['maxOptimizationIterations']], '\n'), con=con, sep='')
  
  # output
  if(is.null(plan[['outputBeamsFile']])) {
    outputBeamsFile <- 'plan.beams'
    plan[['outputBeamsFile']] <- outputBeamsFile
  } else {
    outputBeamsFile <- paste(getwd(), '/', plan[['outputBeamsFile']], sep='')
  }
  
  plan[['voisFile']] <- 'vois.3d'
  plan[['outputValuesFile']] <- 'values.3d'
  
  writeLines(paste('outputBeamsFile = ', outputBeamsFile, '\n'), con=con, sep='')
  writeLines('printSpotPositions = y\n', con=con, sep='')
  writeLines('voisFile = vois.3d\n', con=con, sep='')
  writeLines('outputValuesFile = values.3d\n', con=con, sep='')
  
  close(con)
  
  # chiamata di sistema
  if(outmessages) {ignore.stdout=FALSE; ignore.stderr=FALSE} else {ignore.stdout=TRUE; ignore.stderr=TRUE}
  #dek.setenv <- get('dek.setenv', envir=dektoolsEnv)
  #cmd <- paste('. ', dek.setenv, '; cd ', plan[['name']], '; dek plan.config')
  cmd <- paste('cd ', plan$name, '; dek plan.config')
  system(cmd, ignore.stdout=ignore.stdout, ignore.stderr=ignore.stderr)
  
  save.plan(plan)
  
  return(plan)
}

#' Forward planning (PlanKIT)
#' 
#' Perform a forward planning with PlanKIT.
#' @param plan the Plan object. It should include also the specification of the beams (i.e. the "RTIonPlan").
#' @param outmessages if TRUE, it displays the messages from the PlanKIT output.
#' @return a Plan Object (including data/information of the evaluated output of the forward planning, e.g. values, vois, etc.)
#' 
#' @family PlanKIT
#' @export 
run.dek.forward <- function(plan, outmessages=FALSE) {
  
  message('running plankit forward planning....')
  
  # crea folder
  dir.create(plan[['name']], recursive=TRUE, showWarnings=FALSE)
  
  # rimuove contenuto del folder (il piano viene completamente sovrascritto)
  unlink(paste(plan$name, '/*', sep=''), recursive = FALSE, force = FALSE)
  
  
  # CT e CONTORNI
  if(!is.null(plan[['ct']])) {
    write.3d.array(ct=plan$ct, file.name=paste0(plan$name, '/ct.3d'))
    plan$ctFile <- paste0(plan$name, '/ct.3d')
  }
  if(!is.null(plan[['contours']])) {
    write.contours(contours=plan$contours, name=paste0(plan$name, '/plan'))
    plan$contoursFile <- paste0(plan$name, '/plan.contours')
  }
  
  # BEAMS
  if(!is.null(plan[['beams']])) {
    write.beams(beams=plan[['beams']], format='plankit', file.name=paste0(plan$name, '/plan'))
    plan$inputBeamsFile <- paste0(plan$name, '/plan.beams')
  }
  
  # input files (path assoluti)
  ctFile <- paste(getwd(), '/', plan$ctFile, sep='')
  hounsfieldToDensityFile <- paste(getwd(), '/', plan$hounsfieldToDensityFile, sep='')
  hounsfieldToStoichiometryFile <- paste(getwd(), '/', plan$hounsfieldToStoichiometryFile, sep='')
  contoursFile <- paste(getwd(), '/', plan$contoursFile, sep='')
  plan.config <- paste(plan[['name']], '/plan.config', sep='')
  inputBeamsFile <- paste(getwd(), '/', plan$inputBeamsFile, sep='')
  
  # crea file plan.config
  con <- file(plan.config, "w") # open for writing in text mode
  
  writeLines(paste('ctFile = ', ctFile, '\n'), con=con, sep='')
  writeLines(paste('hounsfieldToDensityFile = ', hounsfieldToDensityFile, '\n'), con=con, sep='')
  writeLines(paste('hounsfieldToStoichiometryFile = ', hounsfieldToStoichiometryFile, '\n'), con=con, sep='')
  writeLines(paste('contoursFile = ', contoursFile, '\n'), con=con, sep='')
  p.computingGridVoxelSizes <- paste(plan[['computingGridVoxelSizes']], collapse=' ')
  writeLines(paste('computingGridVoxelSizes = ', p.computingGridVoxelSizes, '\n', collapse=' '), con=con, sep='')
  writeLines(paste('computingGridCoverage = ', plan[['computingGridCoverage']], '\n'), con=con, sep='')
  writeLines(paste('computingGridCoverageBoundary = ', plan[['computingGridCoverageBoundary']], '\n'), con=con, sep='')
  writeLines(paste('beamLine = ', plan[['beamLine']], '\n'), con=con, sep='')
  writeLines(paste('computingValues = ', plan[['computingValues']], '\n'), con=con, sep='')
  writeLines(paste('cutOffRadius = ', plan[['cutOffRadius']], '\n'), con=con, sep='')
  writeLines(paste('cutEnergyFraction = ', plan[['cutEnergyFraction']], '\n'), con=con, sep='')
  
  # beams input
  if(!is.null(plan$inputBeamsFile)) {
    inputBeamsFile <- paste(getwd(), '/', plan$inputBeamsFile, sep='')
    writeLines(paste('inputBeamsFile = ', inputBeamsFile, '\n'), con=con, sep='')
  }
  
  # beams output
  if(is.null(plan[['outputBeamsFile']])) {
    outputBeamsFile <- 'plan.beams'
    plan[['outputBeamsFile']] <- outputBeamsFile
  } else {
    outputBeamsFile <- paste(getwd(), '/', plan[['outputBeamsFile']], sep='')
  }
  
  plan$voisFile <- 'vois.3d'
  plan$outputValuesFile <- 'values.3d'
  
  writeLines(paste('outputBeamsFile = ', outputBeamsFile, '\n'), con=con, sep='')
  writeLines('printSpotPositions = y\n', con=con, sep='')
  writeLines('voisFile = vois.3d\n', con=con, sep='')
  writeLines('outputValuesFile = values.3d\n', con=con, sep='')
  
  close(con)
  
  # chiamata di sistema
  if(outmessages) {ignore.stdout=FALSE; ignore.stderr=FALSE} else {ignore.stdout=TRUE; ignore.stderr=TRUE}
  #dek.setenv <- get('dek.setenv', envir=dektoolsEnv)
  #cmd <- paste('. ', dek.setenv, '; cd ', plan[['name']], '; dek plan.config')
  cmd <- paste('cd ', plan$name, '; dek plan.config')
  system(cmd, ignore.stdout=ignore.stdout, ignore.stderr=ignore.stderr)
  
  save.plan(plan)
  
  return(plan)
}


# UTILITIES --------------------------------------------------------------------

#' Filepath of objects
#' 
#' Returns the name of the file with complete filepath for the specified object in the plan.
#' @param object the specific object (\code{ctFile}, \code{contoursFile}, \code{hounsfieldToDensityFile}, \code{hounsfieldToStoichiometryFile}, \code{prescriptionFile}, \code{voisFile}, \code{outputValuesFile}, \code{inputBeamsFile}, \code{outputBeamsFile})
#' @return filepah
#' 
#' @family Utilities
get.filepath <- function(object, plan=plan) {
  if(object=='ctFile') {fp <- plan$ctFile}
  else if(object=='contoursFile') {fp <- plan$contoursFile}
  else if(object=='hounsfieldToDensityFile') {fp <- plan$hounsfieldToDensityFile}
  else if(object=='hounsfieldToStoichiometryFile') {fp <- plan$hounsfieldToStoichiometryFile}
  else if(object=='prescriptionFile') {fp <- paste(plan$name, 'plan.prescription', sep='/')}
  else if(object=='voisFile') {fp <- paste(plan$name, 'vois.3d', sep='/')}
  else if(object=='outputValuesFile') {fp <- paste(plan$name, 'values.3d', sep='/')}
  else if(object=='inputBeamsFile') {fp <- plan$inputBeamsFile}
  else if(object=='outputBeamsFile') {fp <- paste(plan$name, plan$outputBeamsFile, sep='/')}
  else{message('no valid file identification.')}
  return(fp)
}
