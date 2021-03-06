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
    'oxygenSimple',
    'protonSimpleGate')
  )
}



#' Available variables (Pure-dek)
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

#' Available constraints (Pure-dek)
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

#' Generate Plan template (Pure-dek)
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
                        hounsfieldToDensityFile=NULL,
                        hounsfieldToStoichiometryFile=NULL,
                        computingGridVoxelSizes=c(2,2,2),
                        computingGridCoverage=0,
                        computingGridCoverageVOI=NULL,
                        computingGridCoverageBoundary=0,
                        beamLines=NULL,
                        computingValues=c('Dose[Gy]'),
                        cutOffRadius=100,
                        cutEnergyFraction=0.01,
                        fields=NULL,
                        prescription=NULL,
                        maxOptimizationIterations=50,
                        outputBeamsFile=NULL,
                        inputBeamsFile=NULL,
                        saveBeamLUTs=FALSE,
                        outputBeamLUTFile=NULL,
                        inputBeamLUTFile=NULL,
                        vois=NULL,
                        outputVoisFile=NULL,
                        inputVoisFile=NULL,
                        values=NULL,
                        outputValuesFile=NULL,
                        inputValuesFile=NULL,
                        beams=NULL)
{

  my.home <- paste0(get.install.path(), '/Pure-dek-install/CT/')
  if(is.null(hounsfieldToDensityFile)) {
    hounsfieldToDensityFile <- paste0(my.home, 'fluka-schneider.HNToDensity.1d')
  }
  if(is.null(hounsfieldToStoichiometryFile)) {
    hounsfieldToStoichiometryFile <- paste0(my.home, 'fluka-schneider.HNToStoichiometry.1d')
  }

  plan <- list(name=name,
               ctFile=ctFile,
               ct=ct,
               contoursFile=contoursFile,
               contours=contours,
               hounsfieldToDensityFile=hounsfieldToDensityFile,
               hounsfieldToStoichiometryFile=hounsfieldToStoichiometryFile,
               computingGridVoxelSizes=computingGridVoxelSizes,
               computingGridCoverage=computingGridCoverage,
               computingGridCoverageVOI=computingGridCoverageVOI,
               computingGridCoverageBoundary=computingGridCoverageBoundary,
               beamLines=beamLines,
               computingValues=computingValues,
               cutOffRadius=cutOffRadius,
               cutEnergyFraction=cutEnergyFraction,
               fields=fields,
               prescription=prescription,
               maxOptimizationIterations=maxOptimizationIterations,
               beams=beams,
               outputBeamsFile=outputBeamsFile,
               inputBeamsFile=inputBeamsFile,
               saveBeamLUTs=saveBeamLUTs,
               outputBeamLUTFile=outputBeamLUTFile,
               inputBeamLUTFile=inputBeamLUTFile,
               vois=vois,
               values=values,
               outputVoisFile=outputVoisFile,
               outputValuesFile=outputValuesFile,
               inputVoisFile=outputVoisFile,
               inputValuesFile=outputValuesFile
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

  message('Saving plan ', plan$name, ' ...')
  plan.file <- paste(plan$name, '/plan.Rdata', sep='')
  dir.create(plan$name, recursive=TRUE, showWarnings=FALSE)
  save(plan, file=plan.file, compress='bzip2')
  cat('plan saved in', plan.file, '\n')

}


#' Read Plan
#'
#' Read a Plan object from disk. The origin folder is derived from the plan name.
#' @param name a string name of the plan. It can be a vector of names. In the latter case a list of plan objects is returned.
#' @return the Plan object (or a list of Plan objects)
#'
#' @family Plan
#' @export
#'
read.plan <- function(name)
{
  if(length(name)==1) {
    message('reading plan: ', name)
    load(paste0(name, '/plan.Rdata'))
    return(plan)
  } else {
    plans <- list()
    for(i in 1:length(name)) {
      message('reading plan: ', name[i])
      load(paste0(name[i], '/plan.Rdata'))
      plans[[i]] <- plan
    }
    return(plans)
  }
}




# PLANKIT ----------------------------------------------------------------------


#' Inverse planning (Pure-dek)
#'
#' Perform an inverse planning with Pure-dek.
#' @param plan the Plan object.
#' @param outmessages if TRUE, it displays the messages from Pure-dek output.
#' @return a Plan Object (including data/information of the evaluated output of the inverse planning, e.g. beams, values, vois, etc.)
#'
#' @family PlanKIT
#' @export
run.dek.inverse <- function(plan, outmessages=FALSE) {

  message('running pure-dek inverse planning...')

  # crea folder
  dir.create(plan$name, recursive=TRUE, showWarnings=FALSE)

  # rimuove contenuto del folder (il piano viene completamente sovrascritto)
  unlink(paste(plan$name, '/*', sep=''), recursive = FALSE, force = FALSE)

  # CT e CONTORNI
  if(!is.null(plan[['ct']])) {
    message('using ct object in plan...')
    write.3d(values=plan[['ct']], file.name=paste0(plan[['name']], '/ct.3d'))
    plan[['ctFile']] <- paste0(plan[['name']], '/ct.3d')
  }
  if(!is.null(plan[['contours']])) {
    write.contours(contours=plan$contours, name=paste0(plan[['name']], '/plan'))
    plan[['contoursFile']] <- paste0(plan[['name']], '/plan.contours')
  }

  # input beamLUT files
  if(!is.null(plan[['inputBeamLUTFile']])) {stop('Usage of external beamLUTs from file not yet implemented.')}

  # BEAMS (carica beams, se sono definiti)
  # assume che l'oggetto beams (se presente) sia il file beams di ingresso.
  if(!is.null(plan[['beams']])) {
    write.beams(beams=plan[['beams']], format='plankit', file.name=paste0(plan[['name']], '/input'))
    plan[['inputBeamsFile']] <- paste0(plan[['name']], '/input.beams')
  } else if(!is.null(plan[['inputBeamsFile']])) {
    plan[['beams']] <- read.beams(plan[['inputBeamsFile']]) # carica comunque i beams e riscrive un file input nuovo
    write.beams(beams=plan[['beams']], format='plankit', file.name=paste0(plan[['name']], '/input'))
    plan[['inputBeamsFile']] <- paste0(plan[['name']], '/input.beams')
  }

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
  plan.config <- paste(plan[['name']], '/plan.config', sep='')
  con <- file(plan.config, "w") # open for writing in text mode

  # CT/Contorni
  writeLines(paste('ctFile = ', plan[['ctFile']], '\n'), con=con, sep='')
  writeLines(paste('hounsfieldToDensityFile = ', plan[['hounsfieldToDensityFile']], '\n'), con=con, sep='')
  writeLines(paste('hounsfieldToStoichiometryFile = ', plan[['hounsfieldToStoichiometryFile']], '\n'), con=con, sep='')
  writeLines(paste('contoursFile = ', plan[['contoursFile']], '\n'), con=con, sep='')

  # calcolo
  p.computingGridVoxelSizes <- paste(plan[['computingGridVoxelSizes']], collapse=' ')
  writeLines(paste('computingGridVoxelSizes = ', p.computingGridVoxelSizes, '\n', collapse=' '), con=con, sep='')
  if(!is.null(plan[['computingGridCoverageVOI']])) {
    coverage <- rep(0, length(plan[['computingGridCoverageVOI']]))
    for(icov in 1:length(coverage)) {
      coverage[icov] <- get.voiindex(voi = plan[['computingGridCoverageVOI']][icov], file.contours=plan[['contoursFile']])
      message('computing grid over ', plan[['computingGridCoverageVOI']][icov], ' (', coverage[icov], ')')
    }
  } else {
    coverage <- paste0(plan[['computingGridCoverage']], collapse=' ')
  }
  writeLines(paste('computingGridCoverage = ', coverage, '\n'), con=con, sep='')
  writeLines(paste('computingGridCoverageBoundary = ', plan[['computingGridCoverageBoundary']], '\n'), con=con, sep='')

  # beamLines (opzionale, ora le beamline sono specificate in field)
  plan[['beamLines']] <- unique(c(as.character(plan[['beamLines']]), as.character(plan[['fields']]$beamLine))) # aggiunge quelle definite in fields
  #if(!is.null(plan[['beamLines']])) {
    # sovrappone direttamente
#     for(indexBL in 1:length(plan[['beamLines']])){
#       message('using beamline: ', plan[['beamLines']][indexBL])
#       writeLines(paste('beamLine = ', plan[['beamLines']][indexBL], '\n'), con=con, sep='')
#     }
  #}

  # computing values
  writeLines('computingValues =', con=con, sep='')
  for(indexCV in 1:length(plan[['computingValues']])) {
    writeLines(paste0(' ', plan[['computingValues']][indexCV]), con=con, sep='')
  }
  writeLines('\n', con=con, sep='')

  writeLines(paste('cutOffRadius = ', plan[['cutOffRadius']], '\n'), con=con, sep='')
  writeLines(paste('cutEnergyFraction = ', plan[['cutEnergyFraction']], '\n'), con=con, sep='')

  # fields
  if(!is.null(plan[['fields']])) { # -> definisce direttamente i fields
    fields <- plan[['fields']]
    if(!is.null(fields)) {
      for(i in 1:nrow(fields)) {
        if(is.na(fields$targetVOIIndex[i])) {
          fields$targetVOIIndex[i] <- get.voiindex(fields$targetVOI[i],
                                                   file.contours=plan[['contoursFile']], contours=plan[['contours']])
        }
        writeLines(paste('beamLine = ', fields[i,'beamLine'], '\n'), con=con, sep='')
        writeLines(paste('targetVOIIndex = ', fields[i,'targetVOIIndex'], '\n'), con=con, sep='')
        writeLines(paste('iecGantryAngle = ', fields[i,'iecGantryAngle'], '\n'), con=con, sep='')
        writeLines(paste('iecPatientSupportAngle = ', fields[i,'iecPatientSupportAngle'], '\n'), con=con, sep='')
        writeLines(paste('interSpotSpacing = ', fields[i,'interSpotSpacing.x'], fields[i,'interSpotSpacing.y'], fields[i,'interSpotSpacing.z'], '\n'), con=con, sep='')
        writeLines(paste('spotsExtensionOutsideTarget = ', fields[i,'spotsExtensionOutsideTarget'], '\n'), con=con, sep='')
        if(!is.na(fields[i,'targetIsocenter.x']) &  !is.na(fields[i,'targetIsocenter.y']) & !is.na(fields[i,'targetIsocenter.z'])) {
          writeLines(paste('isocenter =', fields[i,'targetIsocenter.x'], fields[i,'targetIsocenter.y'], fields[i,'targetIsocenter.z'], '\n'), con=con, sep=' ')
        }
      }
      plan[['fields']] <- fields
    }
  } else { # -> definisce esplicitamente le beamLines da beams
    # beamLines (se si legge da beams, queste devono essere pre-caricate)
    plan[['beamLines']] <- unique(c(as.character(plan[['beamLines']]), as.character(plan[['beams']]$beamLine))) # aggiunge quelle definite in beams
    if(!is.null(plan[['beamLines']])) {
      for(indexBL in 1:length(plan[['beamLines']])){
        writeLines(paste('beamLine = ', plan[['beamLines']][indexBL], '\n'), con=con, sep='')
      }
    }

    # scrive esplicitamente di usare il file degli input beams.
    writeLines(paste('inputBeamsFile = ', plan[['inputBeamsFile']], '\n'), con=con, sep='')
  }

  # optimization
  writeLines(paste('prescriptionFile = ', plan.prescription, '\n'), con=con, sep='')
  writeLines(paste('maxOptimizationIterations = ', plan[['maxOptimizationIterations']], '\n'), con=con, sep='')

  # output BeamLUTs
  if(plan[['saveBeamLUTs']]) {
    plan[['outputBeamLUTFile']] <- paste0(plan[['name']], '/plan')
    writeLines(paste0('outputBeamLUTFile = ', plan[['outputBeamLUTFile']], '\n'), con=con, sep='')
  }

  # output Beams (è fisso)
  plan[['outputBeamsFile']] <- paste0(plan[['name']], '/output.beams')
  writeLines(paste0('outputBeamsFile = ', plan[['outputBeamsFile']], '\n'), con=con, sep='')
  writeLines('printSpotPositions = y\n', con=con, sep='')

  # output Vois (è fisso)
  plan[['outputVoisFile']] <- paste0(plan[['name']], '/vois.3d')
  writeLines(paste0('voisFile = ', plan[['outputVoisFile']], '\n'), con=con, sep='')

  # output Values (è fisso)
  plan[['outputValuesFile']] <- paste0(plan[['name']], '/values.3d')
  writeLines(paste0('outputValuesFile = ', plan[['outputValuesFile']], '\n'), con=con, sep='')

  close(con)

  # chiamata di sistema
  if(outmessages) {ignore.stdout=FALSE; ignore.stderr=FALSE} else {ignore.stdout=TRUE; ignore.stderr=TRUE}
  cmd <- paste0('pure-dek ',  plan[['name']], '/plan.config') # runna da fuori la cartella.
  message('executing: ', cmd)
  t <- system.time(system(cmd, ignore.stdout=ignore.stdout, ignore.stderr=ignore.stderr))
  #message('elapsed time: ', t)
  print(t)
  
  # incolla i beams calcolati nel piano
  plan[['beams']] <- read.beams(beams.file = plan[['outputBeamsFile']])
  save.plan(plan)

  return(plan)
}

#' Forward planning (Pure-dek)
#'
#' Perform a forward planning with Pure-dek. It assumes that the beams are defined in theplan object, checking in this order: beams, inoutBeamsFile, outputBeamsFile.
#' @param plan the Plan object. It should include also the specification of the beams (i.e. the "RTIonPlan").
#' @param outmessages if TRUE, it displays the messages from the Pure-dek output.
#' @return a Plan Object (including data/information of the evaluated output of the forward planning, e.g. values, vois, etc.)
#'
#' @family PlanKIT
#' @export
run.dek.forward <- function(plan, outmessages=FALSE) {

  message('running pure-dek forward planning....')

  # crea folder
  dir.create(plan[['name']], recursive=TRUE, showWarnings=FALSE)

  # rimuove contenuto del folder (il piano viene completamente sovrascritto)
  unlink(paste(plan$name, '/*', sep=''), recursive = FALSE, force = FALSE)

  # CT e CONTORNI
  if(!is.null(plan[['ct']])) {
    message('using ct object in plan...')
    write.3d(values=plan[['ct']], file.name=paste0(plan[['name']], '/ct.3d'))
    plan[['ctFile']] <- paste0(plan[['name']], '/ct.3d')
  }
  if(!is.null(plan[['contours']])) {
    write.contours(contours=plan$contours, name=paste0(plan[['name']], '/plan'))
    plan[['contoursFile']] <- paste0(plan[['name']], '/plan.contours')
  }

  # input beamLUT files
  if(!is.null(plan[['inputBeamLUTFile']])) {stop('Usage of external beamLUTs from file not yet implemented.')}

  # BEAMS (carica beams)
  # assume che l'oggetto beams sia presente
  # Cerca nell'ordine: beams, inputBeamsFile, outputBeamsFile
  if(!is.null(plan[['beams']])) {
    write.beams(beams=plan[['beams']], format='plankit', file.name=paste0(plan[['name']], '/input'))
    plan[['inputBeamsFile']] <- paste0(plan[['name']], '/input.beams')
  } else if(!is.null(plan[['inputBeamsFile']])) {
    plan[['beams']] <- read.beams(plan[['inputBeamsFile']]) # carica comunque i beams e scrive il file input
    write.beams(beams=plan[['beams']], format='plankit', file.name=paste0(plan[['name']], '/input'))
    plan[['inputBeamsFile']] <- paste0(plan[['name']], '/input.beams')
  } else if(!is.null(plan[['outputBeamsFile']])) {
    plan[['beams']] <- read.beams(plan[['outputBeamsFile']]) # carica comunque i beams e scrive il file input
    write.beams(beams=plan[['beams']], format='plankit', file.name=paste0(plan[['name']], '/input'))
    plan[['inputBeamsFile']] <- paste0(plan[['name']], '/input.beams')
  } else {
    stop('beams not defined in plan: ', plan$name)
  }

  # crea file plan.config
  plan.config <- paste(plan[['name']], '/plan.config', sep='')
  con <- file(plan.config, "w") # open for writing in text mode

  # CT e contorni
  writeLines(paste('ctFile = ', plan[['ctFile']], '\n'), con=con, sep='')
  writeLines(paste('hounsfieldToDensityFile = ', plan[['hounsfieldToDensityFile']], '\n'), con=con, sep='')
  writeLines(paste('hounsfieldToStoichiometryFile = ', plan[['hounsfieldToStoichiometryFile']], '\n'), con=con, sep='')
  writeLines(paste('contoursFile = ', plan[['contoursFile']], '\n'), con=con, sep='')

  # calcolo
  p.computingGridVoxelSizes <- paste(plan[['computingGridVoxelSizes']], collapse=' ')
  writeLines(paste('computingGridVoxelSizes = ', p.computingGridVoxelSizes, '\n', collapse=' '), con=con, sep='')
  coverage <- paste0(plan[['computingGridCoverage']], collapse=' ')
  if(!is.null(plan[['computingGridCoverageVOI']])) {
    coverage <- rep(0, length(plan[['computingGridCoverageVOI']]))
    for(icov in 1:length(coverage)) {
      coverage[icov] <- get.voiindex(voi = plan[['computingGridCoverageVOI']][icov], file.contours=plan[['contoursFile']])
      message('computing grid over ', plan[['computingGridCoverageVOI']][icov], ' (', coverage[icov], ')')
    }
  } else {
    coverage <- paste0(plan[['computingGridCoverage']], collapse=' ')
  }
  writeLines(paste('computingGridCoverage = ', coverage, '\n'), con=con, sep='')
  writeLines(paste('computingGridCoverageBoundary = ', plan[['computingGridCoverageBoundary']], '\n'), con=con, sep='')

  # computing values
  writeLines('computingValues =', con=con, sep='')
  for(indexCV in 1:length(plan[['computingValues']])) {
    writeLines(paste0(' ', plan[['computingValues']][indexCV]), con=con, sep='')
  }
  writeLines('\n', con=con, sep='')

  writeLines(paste('cutOffRadius = ', plan[['cutOffRadius']], '\n'), con=con, sep='')
  writeLines(paste('cutEnergyFraction = ', plan[['cutEnergyFraction']], '\n'), con=con, sep='')

  # beamLines (se si legge da beams, queste devono essere pre-caricate)
  plan[['beamLines']] <- unique(c(as.character(plan[['beamLines']]), as.character(plan[['beams']]$beamLine))) # aggiunge quelle definite in beams
  if(!is.null(plan[['beamLines']])) {
    for(indexBL in 1:length(plan[['beamLines']])){
      writeLines(paste('beamLine = ', plan[['beamLines']][indexBL], '\n'), con=con, sep='')
    }
  }

  # beams input
  writeLines(paste('inputBeamsFile = ', plan[['inputBeamsFile']], '\n'), con=con, sep='')

  # output BeamLUTs
  if(plan[['saveBeamLUTs']]) {
    plan[['outputBeamLUTFile']] <- paste0(plan[['name']], '/plan')
    writeLines(paste0('outputBeamLUTFile = ', plan[['outputBeamLUTFile']], '\n'), con=con, sep='')
  }

  # output Beams (è fisso)
  plan[['outputBeamsFile']] <- paste0(plan[['name']], '/output.beams')
  writeLines(paste0('outputBeamsFile = ', plan[['outputBeamsFile']], '\n'), con=con, sep='')
  writeLines('printSpotPositions = y\n', con=con, sep='')

  # output Vois (è fisso)
  plan[['outputVoisFile']] <- paste0(plan[['name']], '/vois.3d')
  writeLines(paste0('voisFile = ', plan[['outputVoisFile']], '\n'), con=con, sep='')

  # output Values (è fisso)
  plan[['outputValuesFile']] <- paste0(plan[['name']], '/values.3d')
  writeLines(paste0('outputValuesFile = ', plan[['outputValuesFile']], '\n'), con=con, sep='')

  close(con)

  # chiamata di sistema
  if(outmessages) {ignore.stdout=FALSE; ignore.stderr=FALSE} else {ignore.stdout=TRUE; ignore.stderr=TRUE}
  cmd <- paste0('pure-dek ',  plan[['name']], '/plan.config') # runna da fuori la cartella.
  message('executing: ', cmd)
  t <- system.time(system(cmd, ignore.stdout=ignore.stdout, ignore.stderr=ignore.stderr))
  #message('elapsed time: ', t)
  print(t)

  # incolla i beams di output (= input) nel piano (già fatto automaticamente)
  # plan[['beams']] <- read.beams(beams.file = plan[['outputBeamsFile']])
  save.plan(plan)

  return(plan)
}


# UTILITIES --------------------------------------------------------------------

#' Embed data into plan
#' 
#' Embed the data (CT, contours, beams, values, vois) associated with a plan directly in the plan object (to provide portability).
#' @param plan A plan.
#' @param embed.ct,embed.contours,embed.beams,embed.values,embed.vois By default all data are embedded.
#' @return A plan with embedded data.
#' @export
#' @family PlanKIT
embed.data <- function(plan, embed.ct=TRUE, embed.contours=TRUE, embed.beams=TRUE, embed.values=TRUE, embed.vois=TRUE)
{
  if(embed.ct & is.null(plan[['ct']])) {plan[['ct']] <- get.ct(plan)}
  if(embed.contours & is.null(plan[['contours']])) {plan[['contours']] <- get.contours(plan)}
  if(embed.beams & is.null(plan[['beams']])) {plan[['beams']] <- get.beams(plan)}
  if(embed.values & is.null(plan[['values']])) {plan[['values']] <- get.values(plan)}
  if(embed.vois & is.null(plan[['vois']])) {plan[['vois']] <- get.vois(plan)}
  return(plan)
}

#' Filepath of objects
#'
#' Returns the name of the file with complete filepath for the specified object in the plan.
#' @param object the specific object (\code{ctFile}, \code{contoursFile}, \code{hounsfieldToDensityFile}, \code{hounsfieldToStoichiometryFile}, \code{prescriptionFile}, \code{voisFile}, \code{outputValuesFile}, \code{inputBeamsFile}, \code{outputBeamsFile})
#' @return filepah
#'
#' @family Utilities
get.filepath <- function(object, plan=plan) {
  if(object=='ctFile') {fp <- plan$ctFile}
  else if(object=='contoursFile') {fp <- plan[['contoursFile']]}
  else if(object=='hounsfieldToDensityFile') {fp <- plan[['hounsfieldToDensityFile']]}
  else if(object=='hounsfieldToStoichiometryFile') {fp <- plan[['hounsfieldToStoichiometryFile']]}
  else if(object=='prescriptionFile') {fp <- paste(plan[['name']], 'plan.prescription', sep='/')}
  else if(object=='voisFile') {fp <- plan[['voisFile']]}
  else if(object=='outputValuesFile') {fp <- plan[['outputValuesFile']]}
  else if(object=='inputBeamsFile') {fp <- plan[['inputBeamsFile']]}
  else if(object=='outputBeamsFile') {fp <- plan[['outputBeamsFile']]}
  else{message('no valid file identification.')}
  return(fp)
}
