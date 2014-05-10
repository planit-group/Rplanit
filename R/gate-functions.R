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
  # da scrivere
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
#' (documentation for the parameters coming soon...)
#' @export
#' @family PlanGate
create.plan.gate <- function(name='gate.simulation',
                             ctFile=NULL,
                             ctWater=NULL,
                             ct=NULL,
                             materialDatabase='./data/MyMaterials.db',
                             HUToMaterialFile='data/ct-HU2mat.txt',
                             contours=NULL,
                             origin=c(0,0,0),    
                             waterIonisationPotential=78,
                             computingValues=c('DoseToMaterial[Gy]'),
                             saveEveryNSeconds=200,
                             size='auto',
                             position=c(0,0,0),
                             computingGridVoxelSizes=c(2,2,2),
                             beams=NULL,
                             beamsFile.gate='./data/beams.gate',
                             particleType='proton',
                             sourceDescriptionFile='./data/ProtonSimple.txt',
                             totalNumberOfPrimaries=0,
                             K=NULL
                             )
{
  plan <- list(name=name,
               ctFile=ctFile,
               ctWater=ctWater, # una semplice lista con la size del waterbox (Dx,Dy,Dz)
               ct=ct,           # possibilità di inserire direttamente un oggetto ct
               materialDatabase=materialDatabase,
               HUToMaterialFile=HUToMaterialFile,
               contours=contours,
               origin=origin, # viene usato solo se si specifica un file *.hdr per la CT. altrimenti si usa direttamente l'origine codificato in ct
               waterIonisationPotential=waterIonisationPotential,
               computingValues=computingValues,
               saveEveryNSeconds=saveEveryNSeconds,
               size=size,
               position=position,
               computingGridVoxelSizes=computingGridVoxelSizes,
               beams=beams, # oggetto beams
               beamsFile.gate=beamsFile.gate, # file dei beams (formato gate)
               particleType=particleType,
               sourceDescriptionFile=sourceDescriptionFile,
               totalNumberOfPrimaries=totalNumberOfPrimaries,
               K=K # pesi di ciascun beam (usato solo per il calcolo dei beam indipendente nelle matrici sparse)    
  )
  class(plan) <- 'gate.plan'
  return(plan)
}


#' Set the CT/phantom (Gate)
#' 
#' Set the patient CT or phantom for a simulation with Gate
#' @param plan.gate a plan (Gate) object.
#' @family PlanGate
set.ct.gate <- function(plan.gate) {
  
  gate.template <- get('gate.template', envir=dektoolsEnv)
  #message('gate.template = ', gate.template)
  
  # legge file Voxel_Patient.mac da template
  vp.mac.template <- file(paste(gate.template, '/mac/Voxel_Patient.mac', sep=''), "rt")
  vp.mac.txt <- readLines(vp.mac.template)
  vp.mac.txt <- paste(vp.mac.txt, collapse='\n')
  close(vp.mac.template)
  
  # OGGETTO CT
  if(!is.null(plan.gate$ct)) {
    message('using ct in plan...')
    write.analyze(values=ct, file.name=paste0(plan.gate$name, '/data/ct'))
    vp.mac.txt <- gsub('#@ct', '', vp.mac.txt)
    vp.mac.txt <- gsub('@myct.hdr', 'data/ct.hdr', vp.mac.txt)
    #plan.gate$origin <- c(ct$x[1], ct$y[1], ct$z[1]) # usa l'origine della CT
    #vp.mac.txt <- sub('@origin', paste(plan.gate$origin, collapse=' '), vp.mac.txt)
  }
  else
  # FILE CT
  if(!is.null(plan.gate$ctFile)) {
    vp.mac.txt <- gsub('#@ct', '', vp.mac.txt)
    vp.mac.txt <- gsub('@myct.hdr', plan.gate$ctFile, vp.mac.txt)
    vp.mac.txt <- sub('@origin', paste(plan.gate$origin, collapse=' '), vp.mac.txt)
  }
  else
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
  #message('gate.template = ', gate.template)

  # crea folder
  dir.create(plan$name, recursive=TRUE, showWarnings=FALSE)
  
  # rimuove le cartelle di gate
  unlink(paste0(plan$name, '/data'), recursive = TRUE, force = FALSE)
  unlink(paste0(plan$name, '/mac'), recursive = TRUE, force = FALSE)
  unlink(paste0(plan$name, '/output'), recursive = TRUE, force = FALSE)

  # crea nuova struttura dal template
  file.copy(paste0(gate.template, '/data'), paste0(plan$name, '/'), recursive=TRUE)
  file.copy(paste0(gate.template, '/mac'), paste0(plan$name, '/'), recursive=TRUE)
  file.copy(paste0(gate.template, '/output'), paste0(plan$name, '/'), recursive=TRUE)
  file.copy(paste0(gate.template, '/run-gate.sh'), paste0(plan$name, '/'))

  # legge file main.mac da template
  main.mac.template <- file(paste(gate.template, '/mac/main.mac', sep=''), "rt")
  main.mac.txt <- readLines(main.mac.template)
  main.mac.txt <- paste(main.mac.txt, collapse='\n')
  close(main.mac.template)
  
  # materiali
  main.mac.txt <- sub('@materialDatabase', plan$materialDatabase, main.mac.txt)
  main.mac.txt <- sub('@waterIonisationPotential', plan$waterIonisationPotential, main.mac.txt)
  
  # main.mac.txt <- sub('@physicsList.gate', plan$physicsList.gate, main.mac.txt)
  main.mac.txt <- gsub('@saveEveryNSeconds', plan$saveEveryNSeconds, main.mac.txt)

  # actor: geometria
  if(plan$size!='auto') {
    main.mac.txt <- sub('#/gate/actor/doseDistribution/setSize', '/gate/actor/doseDistribution/setSize', main.mac.txt)
    main.mac.txt <- sub('@size', paste(plan$size, collapse=' '), main.mac.txt)
  }
  main.mac.txt <- sub('@position', paste(plan$position, collapse=' '), main.mac.txt)
  main.mac.txt <- sub('@voxelSize', paste(plan$computingGridVoxelSizes, collapse=' '), main.mac.txt)

  # actor: dose
  if(sum('DoseToMaterial[Gy]' %in% plan$computingValues)>0) {
    message('computing DoseToMaterial[Gy]')
    main.mac.txt <- sub('@enableDose', 'true', main.mac.txt)
  } else {
    main.mac.txt <- sub('@enableDose', 'false', main.mac.txt)
  }

  if(sum('DoseToMaterial^2[Gy^2]' %in% plan$computingValues)>0) {
    message('computing DoseToMaterial^2[Gy^2]')
    main.mac.txt <- sub('@enableSquaredDose', 'true', main.mac.txt)
  } else {
    main.mac.txt <- sub('@enableSquaredDose', 'false', main.mac.txt)
  }

  if(sum('DoseToMaterialUncertainty[Gy]' %in% plan$computingValues)>0) {
    main.mac.txt <- sub('@enableUncertaintyDose', 'true', main.mac.txt)
  } else {
    main.mac.txt <- sub('@enableUncertaintyDose', 'false', main.mac.txt)
  }

  # da inserire: dose 2 water


  # actor: deposited energy
  if(sum('DepositedEnergy[MeV]' %in% plan$computingValues)>0) {
    main.mac.txt <- sub('@enableEdep', 'true', main.mac.txt)
  } else {
    main.mac.txt <- sub('@enableEdep', 'false', main.mac.txt)
  }

  if(sum('DepositedEnergy^2[MeV^2]' %in% plan$computingValues)>0) {
    main.mac.txt <- sub('@enableSquaredEdep', 'true', main.mac.txt)
  } else {
    main.mac.txt <- sub('@enableSquaredEdep', 'false', main.mac.txt)
  }

  if(sum('DepositedEnergyUncertainty[MeV]' %in% plan$computingValues)>0) {
    main.mac.txt <- sub('@enableUncertaintyEdep', 'true', main.mac.txt)
  } else {
    main.mac.txt <- sub('@enableUncertaintyEdep', 'false', main.mac.txt)
  }

  # actor: number of hits
  if(sum('NumberOfHits' %in% plan$computingValues)>0) {
    main.mac.txt <- sub('@enableNumberOfHits', 'true', main.mac.txt)
  } else {
    main.mac.txt <- sub('@enableNumberOfHits', 'false', main.mac.txt)
  }
  
  # particle type
  main.mac.txt <- sub('@particleType', plan$particleType, main.mac.txt)
  

  # SET BEAMS
  # nota: se non esiste l'oggetto beam occorrerebbe caricarlo e agganciarlo automaticamente...
  # questo sia per la gestione della rotazione del supporto del paziente, sia per il calcolo dell'isocentro...
  if(!is.null(plan$beams)) {
    message('Using beams dataframe stored in plan for ', plan$name)
    file.beams.gate <- paste(plan$name, '/data/beams', sep='')
    write.beams(beams=plan$beams, file.name=file.beams.gate, format='gate')
    
    main.mac.txt <- sub('@beams.gate', './data/beams.gate', main.mac.txt)
    plan$beams.gate <- './data/beams.gate'
  } else {
    main.mac.txt <- sub('@beams.gate', plan$beamsFile.gate, main.mac.txt)
    
    # NON ANCORA IMPLEMENTATO!!!
    # plan$beams <- read.beams.gate(plan)
    stop('beams file reference not yet supported for gate (you have to use directly a beam object: plan$beams <- beams)...')
  }
  
  # check e rotazione supporto del paziente
  patientAngle <- unique(plan$beams$patientAngle)
  if(length(unique(patientAngle))>1) {
    stop('multiple support patient angles not yet implemented in gate...')
  } else {
    main.mac.txt <- sub('@patientSupportAngle', -patientAngle, main.mac.txt)
  }

  # traslazione dell'isocentro
  isocenter <- unique(beams[c('x_iso', 'y_iso', 'z_iso')])
  if(nrow(isocenter)>1) {
    stop('multiple isocenters not yet supported in gate...')
  } else {
    traslazione <- c(sum(range(ct$x))/2-isocenter$x_iso, sum(range(ct$y))/2-isocenter$y_iso, sum(range(ct$z))/2-isocenter$z_iso)
    main.mac.txt <- sub('@translation', paste(traslazione, collapse=' '), main.mac.txt)
  }
  
  # numero di eventi
  #main.mac.txt <- sub('@totalNumberOfPrimaries', plan$totalNumberOfPrimaries, main.mac.txt)
  
  # scrive file
  #print(main.mac.txt)
  main.mac <- file(paste(plan$name, '/mac/main.mac', sep=''), "wt")
  write(main.mac.txt, file=main.mac)
  close(main.mac)
  
  # SET CT
  set.ct.gate(plan)
  
  # SET RUN
  run.mac <- file(paste(plan$name, '/mac/run.mac', sep=''), "wt")
  run.mac.txt <- paste('/gate/application/setTotalNumberOfPrimaries', plan$totalNumberOfPrimaries)
  write(run.mac.txt, file=run.mac)
  close(run.mac)
  
  message('gate structure ', plan$name, ' written.')
  return(plan)
}


#' Evaluate forward planning (Gate)
#' 
#' @param plan The plan object (Gate)
#' @param N The mumber of primary particles to simulate (events). If specified it overwrites the number specified in the plan.
#' @param evaluate.sparse.array Evaluate and save also a sparse array for the values (in which the information for the individual beam contribution is stored).
#' @param outmessages Show the optuputs from the Gate simulation
#' @return The updated plan object (Gate)
#' @export
#' @family PlanGate
run.gate.forward <- function(plan=plan, N=NULL, K=NULL, evaluate.sparse.arrays=FALSE, outmessages=FALSE)
{
  
  #gate.template <- get('gate.template', envir=dektoolsEnv)
  #gate.setenv <- get('gate.setenv', envir=dektoolsEnv)
  
  if(!is.null(N)) {
    message('using N=', N, ' primary particles...')
    plan$totalNumberOfPrimaries <- N
  }
  if(!is.null(K)) {
    message('using K...')
    plan$K <- K
  }
  
  # crea struttura file
  plan <- create.gate.structure(plan)
  
  # SET bash script
  #run.template <- file(paste(gate.template, '/run-gate.sh', sep=''), "rt")
  #run.txt <- readLines(run.template)
  #run.txt <- paste(run.txt, collapse='\n')
  #close(run.template)
  #run.txt <- sub('@gate.setenv', gate.setenv, run.txt)
  #run <- file(paste(plan$name, '/run-gate.sh', sep=''), "wt")
  #write(run.txt, file=run)
  #close(run)
  
  # run
  cmd <- paste('cd ', plan$name, '; chmod +x ./run-gate.sh; ./run-gate.sh', sep='')
  if(outmessages) {ignore.stdout=FALSE; ignore.stderr=FALSE} else {ignore.stdout=TRUE; ignore.stderr=TRUE}
  
  if(!evaluate.sparse.arrays) {
    message('running Gate...')
    system(cmd, ignore.stdout=ignore.stdout, ignore.stderr=ignore.stderr)
  }
  else 
  {
    #recupera beams
    beams <- get.beams(plan)
    beams.file <- plan$beamsFile.gate # il file deve essere comunque presente (da create.gate.structure)
    if(beams.file=='./data/beams.gate') {beams.file <- paste0(plan$name, '/data/beams.gate')}
    Nb <- nrow(beams)
    N <- plan$totalNumberOfPrimaries
    if(is.null(plan$K)) {message('using default K = cost. ...'); plan$K <- rep(1, Nb)}
    f <- sum(plan$K*beams$fluence)/N
    plan$K <- plan$K/f
    Ne <- round(beams$fluence*plan$K)
    #print(Ne)
    
    for(b in 1:Nb){
      message('evaluating beam ', b, '/', Nb, ' (', format(b/Nb*100, digits=3), '%) --- ',
              'N primaries: ', Ne[b], '/', N , ' (', format(Ne[b]/N*100, digits=3), '%) ...')
      
      if(Ne[b]==0) next
      # beams
      beam <- beams[b,]
      write.beams(beams=beam, file.name=beams.file, format='gate', add.extension=FALSE)
      
      # numero di eventi
      run.mac <- file(paste(plan$name, '/mac/run.mac', sep=''), "wt")
      run.mac.txt <- paste('/gate/application/setTotalNumberOfPrimaries', Ne[b])
      write(run.mac.txt, file=run.mac)
      close(run.mac)
      
      # run
      system(cmd, ignore.stdout=ignore.stdout, ignore.stderr=ignore.stderr)
      
      vb <- get.values.gate.plan(plan)
      vs.tmp <- get.sparse.array.from.values(vb)
      vs.tmp$beam.id <- b
      vs.tmp$Ne <- Ne[b]
      
      if(b==1) {values.sparse <- vs.tmp} else {values.sparse <- rbind(values.sparse, vs.tmp)}
    }
    
    message('saving sparse matrix...')
    save(values.sparse, file=paste(plan$name, 'values.sparse.Rdata', sep='/'), compress='bzip2')
  }
  message('...done.')
  
  # salva oggetto plan
  save.plan(plan)
  
  return(plan)
}
