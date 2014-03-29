#' Available variables (Gate)
#' 
#' The available variables name identifiers that can be evaluated with Gate.
#' @return A list of the available variables
#' @export
available.variables.gate <- function()
{
  return(
    c('Dose[Gy]',
      'Dose^2[Gy^2]',
      'DoseUncertainty[Gy]',
      'DoseToMaterial[Gy]',
      'DoseToMaterial^2[Gy^2]',
      'DoseToMaterialUncertainty[Gy]',
      'DepositedEnergy[MeV]',
      'DepositedEnergy^2[MeV^2]',
      'DepositedEnergyUncertainty[MeV]',
      'NumberOfHits')
    )
}


#' Map the variable name to file name (Gate)
#' @param variable The variable name.
#' @return The prefix of the file name.
#' @export
#' @family PlanGate
file.variable.gate <- function(variable)
{
  file.names <- c('values-Dose', # dose 2 water
                  'values-Dose-Squared',
                  'values-Dose-Uncertainty',
                  'values-Dose', # dose 2 material
                  'values-Dose-Squared',
                  'values-Dose-Uncertainty',
                  'values-Edep',
                  'values-Edep-Squared',
                  'values-Dose-Uncertainty',
                  'values-NbOfHits')
  variables <- available.variables.gate()
  return(file.names[which(variables==variable)])
}

# GATE MATERIALS ---------------------------------------------------------------

#' Read materials database (Gate)
#' 
#' @param materialDatabase the file name of the material database
#' @return a list of materials
#' @export
#' @family GateMaterials
read.materials.gate <- function(materialDatabase=NULL)
{
  if(is.null(materialDatabase)) {
    # read directly from the standard Gate material database
  }
}

#' Generate Materials (Gate)
#' 
#' Generate the material database for Gate simulations.
#' @param MaterialTable name of the file containing the conversion of CT number to elemental weight.
#' @param DensityTable name of the file containing the conversion of CT number to density values.
#' @param DensityTolerance the density step for the generation of the material database.
#' @param OutputMaterialDatabase output file name for the material database.
#' @param OutputHUMaterial output file name for the CT number to material conversion.
#' @export
#' @family GateMaterials
generate.materials.gate <- function(MaterialTable='./data/Schneider2000MaterialsTable.txt',
                                    DensityTable='./data/Schneider2000DensitiesTable.txt',
                                    DensityTolerance=0.1,
                                    OutputMaterialDatabase='./data/ct-HUmaterials.db',
                                    OutputHUMaterial='./data/ct-HU2mat.txt')
{
  
}


#' Create a water box object (Gate)
#' 
#' Create a water box object to be used as a CT with Gate simulations.
#' @param Dx The size of the box along x (mm).
#' @param Dy The size of the box along y (mm).
#' @param Dz The size of the box along z (mm).
#' @export
#' @family PlanGate
create.waterbox.gate <- function(Dx=100, Dy=100, Dz=100)
{
  return(list(Dx=Dx, Dy=Dy, Dz=Dz))
}

# PLAN GATE --------------------------------------------------------------------

#' Generate Plan template (Gate)
#' 
#' Generate a default (empty) template of the Plan, to be used with Gate.
#' @export
#' @family PlanGate
create.plan.gate <- function(name='gate.simulation',
                             ctFile.gate=NULL,
                             ctWater=NULL,
                             materialDatabase.gate='./data/MyMaterials.db',
                             HUToMaterialFile='data/ct-HU2mat.txt',
                             origin.gate=c(0,0,0),    
                             waterIonisationPotential.gate=78,
			                       computingValues.gate='DoseToMaterial[Gy]',
                             saveEveryNSeconds.gate=200,
                             size.gate='auto',
                             position.gate=c(0,0,0),
                             voxelSize.gate=c(2,2,2),
			                       beams=NULL,
                             beams.gate='./data/beams.gate',
                             particleType.gate='proton',
                             sourceDescriptionFile.gate='./data/ProtonSimple.txt',
                             totalNumberOfPrimaries=0
                             )
{
  plan <- list(name=name,
               ctFile.gate=ctFile.gate,
               ctWater=ctWater, # una semplice lista con la size del waterbox (Dx,Dy,Dz)
               materialDatabase.gate=materialDatabase.gate,
               HUToMaterialFile=HUToMaterialFile,
               origin.gate=origin.gate,
               waterIonisationPotential.gate=waterIonisationPotential.gate,
	             computingValues.gate=computingValues.gate,
               saveEveryNSeconds.gate=saveEveryNSeconds.gate,
               size.gate=size.gate,
               position.gate=position.gate,
               voxelSize.gate=voxelSize.gate,
	             beams=beams, # oggetto beams
               beams.gate=beams.gate, # file dei beams
               particleType.gate=particleType.gate,
               sourceDescriptionFile.gate=sourceDescriptionFile.gate,
               totalNumberOfPrimaries=totalNumberOfPrimaries    
  )
  return(plan)
}


#' Set the CT/phantom (Gate)
#' 
#' Set the patient CT or phantom for a simulation with Gate
#' @param plan.gate a plan (Gate) object.
#' @family PlanGate
set.ct.gate <- function(plan.gate) {
  
  gate.template <- get('gate.template', envir=dektoolsEnv)
  message('gate.template = ', gate.template)
  
  # legge file Voxel_Patient.mac da template
  vp.mac.template <- file(paste(gate.template, '/mac/Voxel_Patient.mac', sep=''), "rt")
  vp.mac.txt <- readLines(vp.mac.template)
  vp.mac.txt <- paste(vp.mac.txt, collapse='\n')
  close(vp.mac.template)
  
  # WATERBOX
  if(!is.null(plan.gate$ctWater)) {
    message('setting water box...')
    vp.mac.txt <- gsub('#@wb', '', vp.mac.txt)
    vp.mac.txt <- sub('@Dx', plan.gate$ctWater$Dx, vp.mac.txt)
    vp.mac.txt <- sub('@Dy', plan.gate$ctWater$Dy, vp.mac.txt)
    vp.mac.txt <- sub('@Dz', plan.gate$ctWater$Dz, vp.mac.txt)
  } else {
    error('setting for CT not yet implemented')
  }
  
  # write file .mac
  vp.mac <- file(paste(plan.gate$name, '/mac/Voxel_Patient.mac', sep=''), "wt")
  write(vp.mac.txt, file=vp.mac)
  close(vp.mac)
  message('CT/phantom written in plan ', plan.gate$name)
}


#' Create Gate structure
#' 
#' Create the folder structure of a Gate simulation for a specific plan.
#' @param plan A Plan Object.
#' @export
#' @family PlanGate
create.gate.structure <- function(plan)
{
  gate.template <- get('gate.template', envir=dektoolsEnv)
  message('gate.template = ', gate.template)

  # crea folder
  dir.create(plan$name, recursive=TRUE, showWarnings=FALSE)
  
  # rimuove le cartelle di gate
  unlink(paste(plan$name, '/data', sep=''), recursive = TRUE, force = FALSE)
  unlink(paste(plan$name, '/mac', sep=''), recursive = TRUE, force = FALSE)
  unlink(paste(plan$name, '/output', sep=''), recursive = TRUE, force = FALSE)

  # crea nuova struttura dal template
  file.copy(paste(gate.template, '/data', sep=''), paste(plan$name, '/', sep=''), recursive=TRUE)
  file.copy(paste(gate.template, '/mac', sep=''), paste(plan$name, '/', sep=''), recursive=TRUE)
  file.copy(paste(gate.template, '/output', sep=''), paste(plan$name, '/', sep=''), recursive=TRUE)

  # legge file main.mac da template
  main.mac.template <- file(paste(gate.template, '/mac/main.mac', sep=''), "rt")
  main.mac.txt <- readLines(main.mac.template)
  main.mac.txt <- paste(main.mac.txt, collapse='\n')
  close(main.mac.template)
  
  # materiali
  main.mac.txt <- sub('@materialDatabase.gate', plan$materialDatabase.gate, main.mac.txt)
  main.mac.txt <- sub('@waterIonisationPotential.gate', plan$waterIonisationPotential.gate, main.mac.txt)
  
  # main.mac.txt <- sub('@physicsList.gate', plan$physicsList.gate, main.mac.txt)
  main.mac.txt <- gsub('@saveEveryNSeconds.gate', plan$saveEveryNSeconds.gate, main.mac.txt)

  # actor: geometria
  if(plan$size.gate!='auto') {
    main.mac.txt <- sub('#/gate/actor/doseDistribution/setSize', '/gate/actor/doseDistribution/setSize', main.mac.txt)
    main.mac.txt <- sub('@size.gate', paste(plan$size.gate, collapse=' '), main.mac.txt)
  }
  main.mac.txt <- sub('@position.gate', paste(plan$position.gate, collapse=' '), main.mac.txt)
  main.mac.txt <- sub('@voxelSize.gate', paste(plan$voxelSize.gate, collapse=' '), main.mac.txt)

  # actor: dose
  if(sum('DoseToMaterial[Gy]' %in% plan$computingValues.gate)>0) {
    main.mac.txt <- sub('@enableDose', 'true', main.mac.txt)
  } else {
    main.mac.txt <- sub('@enableDose', 'false', main.mac.txt)
  }

  if(sum('DoseToMaterial^2[Gy^2]' %in% plan$computingValues.gate)>0) {
    main.mac.txt <- sub('@enableSquaredDose', 'true', main.mac.txt)
  } else {
    main.mac.txt <- sub('@enableSquaredDose', 'false', main.mac.txt)
  }

  if(sum('DoseToMaterialUncertainty[Gy]' %in% plan$computingValues.gate)>0) {
    main.mac.txt <- sub('@enableUncertaintyDose', 'true', main.mac.txt)
  } else {
    main.mac.txt <- sub('@enableUncertaintyDose', 'false', main.mac.txt)
  }

  # da inserire: dose 2 water


  # actor: deposited energy
  if(sum('DepositedEnergy[MeV]' %in% plan$computingValues.gate)>0) {
    main.mac.txt <- sub('@enableEdep', 'true', main.mac.txt)
  } else {
    main.mac.txt <- sub('@enableEdep', 'false', main.mac.txt)
  }

  if(sum('DepositedEnergy^2[MeV^2]' %in% plan$computingValues.gate)>0) {
    main.mac.txt <- sub('@enableSquaredEdep', 'true', main.mac.txt)
  } else {
    main.mac.txt <- sub('@enableSquaredEdep', 'false', main.mac.txt)
  }

  if(sum('DepositedEnergyUncertainty[MeV]' %in% plan$computingValues.gate)>0) {
    main.mac.txt <- sub('@enableUncertaintyEdep', 'true', main.mac.txt)
  } else {
    main.mac.txt <- sub('@enableUncertaintyEdep', 'false', main.mac.txt)
  }

  # actor: number of hits
  if(sum('NumberOfHits' %in% plan$computingValues.gate)>0) {
    main.mac.txt <- sub('@enableNumberOfHits', 'true', main.mac.txt)
  } else {
    main.mac.txt <- sub('@enableNumberOfHits', 'false', main.mac.txt)
  }

  # beams
  if(!is.null(plan$beams)) {
    message('Using beams dataframe stored in plan for ', plan$name)
    file.beams.gate <- paste(plan$name, '/data/beams', sep='')
    write.beams(beams=plan$beams, file.name=file.beams.gate, format='gate')
    main.mac.txt <- sub('@beams.gate', './data/beams.gate', main.mac.txt)
    plan.gate$beams.gate <- './data/beams.gate'
  } else {
    main.mac.txt <- sub('@beams.gate', plan$beams.gate, main.mac.txt)
  }
  
  # numero di eventi
  main.mac.txt <- sub('@totalNumberOfPrimaries', plan$totalNumberOfPrimaries, main.mac.txt)
  
  # scrive file
  #print(main.mac.txt)
  main.mac <- file(paste(plan$name, '/mac/main.mac', sep=''), "wt")
  write(main.mac.txt, file=main.mac)
  close(main.mac)
  message('plan ', plan$name, ' written.')
  
  # SET CT
  set.ct.gate(plan)
  
  return(plan)
}


#' calcola una simulazione Forward Planning
run.gate.forward <- function(plan=plan, N=plan$TotalNumberOfPrimaries, save.sparse.arrays=FALSE, outmessages=FALSE)
{
  
  gate.template <- get('gate.template', envir=dektoolsEnv)
  gate.setenv <- get('gate.setenv', envir=dektoolsEnv)
  
  plan$totalNumberOfPrimaries <- N
  
  # crea struttura file
  plan <- create.gate.structure(plan)
  
  # SET bash script
  run.template <- file(paste(gate.template, '/run-gate.sh', sep=''), "rt")
  run.txt <- readLines(run.template)
  run.txt <- paste(run.txt, collapse='\n')
  close(run.template)
  run.txt <- sub('@gate.setenv', gate.setenv, run.txt)
  run <- file(paste(plan$name, '/run-gate.sh', sep=''), "wt")
  write(run.txt, file=run)
  close(run)
  
  # run
  message('running Gate...')
  if(outmessages) {ignore.stdout=FALSE; ignore.stderr=FALSE} else {ignore.stdout=TRUE; ignore.stderr=TRUE}
  cmd <- paste('cd ', plan$name, '; chmod +x ./run-gate.sh; ./run-gate.sh', sep='')
  system(cmd, ignore.stdout=ignore.stdout, ignore.stderr=ignore.stderr)
  message('...done.')
  
  return(plan)
}
