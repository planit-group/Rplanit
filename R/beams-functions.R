# BEAMS ------------------------------------------------------------------------

#' Crea oggetto field (dataframe ops)
#' 
#' nota: si tratta di un dataframe: diversi fields sono semplicemente aggiunti
#' appendendo delle righe al dataframe...
#' 
#' @family Beams
#' @export
#' 
create.field <- function(targetVOI='PTV',
                         targetVOIIndex=NA,
                         iecGantryAngle=0,
                         iecPatientSupportAngle=0,
                         interSpotSpacing=c(1.5,1.5,1),
                         spotsExtensionOutsideTarget=0,
                         targetIsocenter.x=NA,
                         targetIsocenter.y=NA,
                         targetIsocenter.z=NA) 
{  
  field <- data.frame(targetVOI=targetVOI,
                      targetVOIIndex=targetVOIIndex,
                      iecGantryAngle=iecGantryAngle,
                      iecPatientSupportAngle=iecPatientSupportAngle,
                      interSpotSpacing.x=interSpotSpacing[1],
                      interSpotSpacing.y=interSpotSpacing[2],
                      interSpotSpacing.z=interSpotSpacing[3],
                      spotsExtensionOutsideTarget=spotsExtensionOutsideTarget,
                      targetIsocenter.x=targetIsocenter.x,
                      targetIsocenter.y=targetIsocenter.y,
                      targetIsocenter.z=targetIsocenter.z)
  return(field)
}


#' Read beams from file (PlanKIT)
#' 
#' Read RtIonPlan file (Plankit *.beams format) and return a beams object.
#' @param beams.file The file name.
#' @return a beams object.
#' @family Beams
#' @export
read.beams <- function(beams.file)
{
  message('reading beams file...')
  beams <- read.table(beams.file, skip=1)
  
  # mette i nomi delle colonne
  if(ncol(beams)==12) {
    names(beams) <- c('x_iso', 'y_iso', 'z_iso', 'gantryAngle', 'patientAngle', 'fluence', 'energy', 'deflX', 'deflY', 'x_s', 'y_s', 'z_s')
    message('Spot positions present...')
  } else if(ncol(beams)==9) {
    names(beams) <- c('x_iso', 'y_iso', 'z_iso', 'gantryAngle', 'patientAngle', 'fluence', 'energy', 'deflX', 'deflY')
    message('Spot positions NOT present...')
  }
  
  return(beams)
}

#' Read beams from file (Fluka)
#' 
#' Read RtIonPlan file (Fluka format) and return a beams object.
#' @param fluka.file The file name.
#' @return a beams object.
#' @family Beams
#' @export
read.beams.fluka <- function(fluka.file)
{
  message('reading beams file...')
  beams.fluka <- readLines(fluka.file)
  index <- 0
  Nbeams <- 0
  
  # ciclo preliminare per numero totale di fasci
  for(i in 1:length(beams.fluka)) {
    if(length(grep('#points', beams.fluka[i], fixed=TRUE))>0) {
      st <- strsplit(beams.fluka[i], ' +')[[1]]
      Nbeams <- Nbeams + as.numeric(st[2])
    }
  }
  message('number of beams: ', Nbeams)
  
  beams <- create.beam(nbeams=Nbeams, with.spots=FALSE)
  
  i <- 0
  while(i < length(beams.fluka)) {
    i <- i+1
    if(length(grep('submachine#', beams.fluka[i], fixed=TRUE))>0) {
      st <- strsplit(beams.fluka[i], ' +')[[1]]
      energy <- as.numeric(st[3])
      message(i, ', energy: ', energy)
    } else if (length(grep('#points', beams.fluka[i], fixed=TRUE))>0) {
      while(
        length(grep('submachine#', beams.fluka[i], fixed=TRUE))==0 & i < length(beams.fluka)
        ) {
        i <- i+1
        if(length(grep('submachine#', beams.fluka[i], fixed=TRUE))==0) {
          st <- strsplit(beams.fluka[i], ' +')[[1]]
          deflX <- as.numeric(st[1])
          deflY <- as.numeric(st[2])
          fluence <- as.numeric(st[3])
          index <- index+1
          beams$energy[index] <- energy
          beams$deflX[index] <- deflX
          beams$deflY[index] <- deflY
          beams$fluence[index] <- fluence
          if(is.na(deflX)) {message('NA: ', i)}
        }
        #print(beams[index,])
      }
      i <- i-1
    }
  }

  return(beams)
}


#' recupera oggetto beams
#' 
#' input serve a selezionare i beams usati come "input" per il calcolo
#' di default si recuperano quelle in "output".
#' in teoria sarebbero identici tranne casi particolari per il posizionamento
#' degli spot (se è stata cambiata CT tra inv e fwd planning).
#' 
#' @family Beams
#' @export
get.beams <- function(plan, input=FALSE)
{
  #if(!is.null(plan$outputBeamsFile)) {
  #  beams.file <- get.filepath('outputBeamsFile', plan=plan)
  #} else if(!is.null(plan$inputBeamsFile)) {
  #  beams.file <- get.filepath('inputBeamsFile', plan=plan)
  #}
  
  if(input){beams.file <- get.filepath('inputBeamsFile', plan=plan)}
  else {beams.file <- get.filepath('outputBeamsFile', plan=plan)}
  return(read.beams(beams.file))
}


#' converte fluenze in beams da MU a numero di particelle
#' 
#' secondo le curve di calibrazione specificate
#' 
#' @family Beams
#' @export
convert.mu2fluence <- function(beams, type='protonCnao')
{
  if(type=='protonCnao') {
    E <- c(62.28, 80.98, 96.84, 118.19, 147.72, 173.61, 196.68, 218.55, 226.91)
    m2f <- c(0.9656, 0.9581, 0.9555, 0.9680, 0.9731, 0.9760, 0.9788, 0.9886, 0.9900)
  } else if(type=='carbonCnao') {
    E <- c(115.23, 150.71, 181.17, 208.58, 257.5, 279.97, 322.08, 332.15, 380.45, 398.84)
    m2f <- c(0.945, 0.9526, 0.9566, 0.9598, 0.9613, 0.9645, 0.9626, 0.9661, 0.9667, 0.97)
  }
  
  mu2fluence <- approxfun(x=E, y=m2f)
  beams$fluence <- beams$fluence * mu2fluence(beams$energy)
  return(beams)
}


#' converte fluenze in beams da numero di particelle a MU
#' 
#' secondo le curve di calibrazione specificate
#' 
#' @family Beams
#' @export
convert.fluence2mu <- function(beams, type='protonCnao')
{
  if(type=='protonCnao') {
    E <- c(62.28, 80.98, 96.84, 118.19, 147.72, 173.61, 196.68, 218.55, 226.91)
    f2m <- c(0.9656, 0.9581, 0.9555, 0.9680, 0.9731, 0.9760, 0.9788, 0.9886, 0.9900)^-1
  } else if(type=='carbonCnao') {
    E <- c(115.23, 150.71, 181.17, 208.58, 257.5, 279.97, 322.08, 332.15, 380.45, 398.84)
    f2m <- c(0.945, 0.9526, 0.9566, 0.9598, 0.9613, 0.9645, 0.9626, 0.9661, 0.9667, 0.97)^-1
  }
  
  mu2fluence <- approxfun(x=E, y=f2m)
  beams$fluence <- beams$fluence * mu2fluence(beams$energy)
  return(beams)
}


#' Create an empty beams object
#' 
#' @param nbeams The number of the beams
#' @param with.spots Include spot coordinates
#' @family Beams
#' @export
create.beam <- function(nbeams=1, with.spots=FALSE)
{
  if(with.spots) {
    beam <- data.frame(x_iso=rep(0, nbeams), y_iso=0, z_iso=0, gantryAngle=0, patientAngle=0, fluence=1, energy=0, deflX=0, deflY=0, x_s=0, y_s=0, z_s=0)
  } else {
    beam <- data.frame(x_iso=rep(0, nbeams), y_iso=0, z_iso=0, gantryAngle=0, patientAngle=0, fluence=1, energy=0, deflX=0, deflY=0)
  }
  return(beam)
}


#' Write beams (RTIonPlan)
#' 
#' Write the beams object on file. It can use different formats (default PlanKIT).
#' 
#' @param beams The beams object
#' @param file.name The prefix of the file name.
#' @param format The format. It can be 'plankit', 'gate' or 'fluka'.
#' @param ion. The primary ion. (used by fluka).
#' @family Beams
#' @export
#' 
write.beams <- function(beams, file.name, format='plankit', ion='1H')
{
  
  # PlanKIT
  if(format=='plankit') {
    file.name <- paste(file.name, '.beams', sep='')
    N <- nrow(beams)
    con <- file(file.name, "w")
    writeLines(paste(N), con=con)
    close(con)
    write.table(beams, file=file.name, append=TRUE, col.names=FALSE, row.names=FALSE)
  }
  
  # Gate
  if(format=='gate') {
    
    # variabili standard non contenute in *.beams
    PlanName <- "DEK-IMPT"
    NumberOfFraction <- 1
    FractionID <- 1
    SpotTunnedID <- 4 # che d'è???
    
    # evaluate number of fields from unique gantry and patient angles
    angles <- unique(cbind(beams$gantryAngle, beams$patientAngle))
    GantryAngle <- angles[,1]
    PatientSupportAngle <- angles[,2]
    NumberOfFields <- dim(angles)[1]
    
    # somma dei pesi/fluenze/MU
    TotalMetersetWeightOfAllFields <- sum(beams$fluence)
    FinalCumulativeMeterSetWeight <- aggregate(beams$fluence, by=list(beams$gantryAngle, beams$patientAngle), FUN='sum')$x # MU erogate in totale (MU del fascio)
    CumulativeMeterSetWeight <- 0 # MU erogate a quel punto
    
    # posizione isocentro (assume unica per ogni fields)
    IsocenterPosition <- aggregate(cbind(beams$x_iso, beams$y_iso, beams$z_iso), by=list(beams$gantryAngle, beams$patientAngle), FUN='mean')[,3:5]
    
    # energie per field
    NumberOfEnergies <- rep(0, NumberOfFields)
    for (i in 1:NumberOfFields) {
      NumberOfEnergies[i] <- length(unique(beams$energy[beams$gantryAngle==GantryAngle[i] & beams$patientAngle==PatientSupportAngle[i]]))
    }
    NumberOfControlPoints <- NumberOfEnergies + 1 # forse si assume che il secondo control point di una slice sia coincidente con il primo della slice successiva...
    
    
    # scrivi file *.gate
    
    # open connection to file
    gate <- file(paste(file.name, '.gate', sep=''), "w")
    
    cat('#TREATMENT-PLAN-DESCRIPTION\n', file=gate, sep='')
    cat('#PlanName\n', PlanName, '\n', file=gate, sep='')
    
    cat('#NumberOfFractions\n', NumberOfFraction, '\n',file=gate, sep='')
    cat('#FractionID\n', FractionID, '\n',file=gate, sep='')
    
    cat('#NumberOfFields\n', NumberOfFields, '\n',file=gate, sep='')
    for (i in 1:NumberOfFields) {
      cat('#FieldsID\n', i, '\n',file=gate, sep='')
    }
    cat('#TotalMetersetWeightOfAllFields\n', TotalMetersetWeightOfAllFields, '\n',file=gate, sep='')
    
    # fields main loop
    for (i in 1:NumberOfFields) {
      cat('\n#FIELD-DESCRIPTION\n', file=gate, sep='')
      cat('#FieldID\n', i, '\n',file=gate, sep='')
      cat('#FinalCumulativeMeterSetWeight\n', FinalCumulativeMeterSetWeight, '\n',file=gate, sep='')
      cat('#GantryAngle\n', GantryAngle[i], '\n',file=gate, sep='')
      cat('#PatientSupportAngle\n', PatientSupportAngle[i], '\n',file=gate, sep='')
      cat('#IsocenterPosition\n', file=gate)
      cat(IsocenterPosition[i,1], IsocenterPosition[i,2], IsocenterPosition[i,3], '\n',file=gate)
      cat('#NumberOfControlPoints\n', NumberOfControlPoints[i], '\n',file=gate, sep='')
      
      cat('\n#SPOTS-DESCRIPTION\n', file=gate, sep='')
      
      CumulativeMeterSetWeight <- 0
      Energy <- unique(beams$energy[beams$gantryAngle==GantryAngle[i] & beams$patientAngle==PatientSupportAngle[i]])
      
      # ciclo sulle energie/punti di controllo
      for (j in 1:(NumberOfEnergies[i])) { 
        cat('#ControlPointIndex\n', j-1, '\n',file=gate, sep='')
        cat('#SpotTunnedID\n', SpotTunnedID, '\n',file=gate, sep='')
        cat('#CumulativeMeterSetWeight\n', CumulativeMeterSetWeight, '\n',file=gate, sep='')
        cat('#Energy (MeV)\n', Energy[j], '\n',file=gate, sep='')    
        Spots <- beams[beams$gantryAngle==GantryAngle[i] & beams$patientAngle==PatientSupportAngle[i] & beams$energy==Energy[j],]
        NbOfScannedSpots <- dim(Spots)[1]
        cat('#NbOfScannedSpots\n', NbOfScannedSpots, '\n',file=gate, sep='')
        cat('#X Y Weight\n',file=gate, sep='')
        for (k in 1:NbOfScannedSpots) {
          cat(Spots$deflX[k], Spots$deflY[k], Spots$fluence[k], '\n', file=gate)
        }    
        # aggiunge il contributo della slice alle MU/# particelle totali.
        CumulativeMeterSetWeight <- CumulativeMeterSetWeight + sum(Spots$fluence)    
      }
      
      # nel formatto c'è un'ultima slice "dummy" con energia=0 e NbOfScannedSpots=0
      cat('#ControlPointIndex\n', j, '\n',file=gate, sep='')
      cat('#SpotTunnedID\n', 0, '\n',file=gate, sep='')
      cat('#CumulativeMeterSetWeight\n', CumulativeMeterSetWeight, '\n',file=gate, sep='')
      cat('#Energy (MeV)\n', 0, '\n',file=gate, sep='')       
      cat('#NbOfScannedSpots\n', 0, '\n',file=gate, sep='')
      cat('#X Y Weight\n',file=gate, sep='')
    }
    
    close(gate)
  }
  
  
  # Fluka
  if(format=='fluka') {
    
    # header
    ions <- c('1H', '12C')
    charges <- c(1, 6)
    masses <- c(1, 12)
    bolus <- 0
    ripplefilter <- 0
    patient_id <- 1
    charge <- charges[which(ions==ion)]
    mass <- masses[which(ions==ion)]
    stepsize <- '1.0 1.0'
    
    # evaluate number of fields from unique gantry and patient angles
    angles <- unique(cbind(beams$gantryAngle, beams$patientAngle)) 
    GantryAngle <- angles[,1]
    PatientSupportAngle <- angles[,2]
    NumberOfFields <- dim(angles)[1]
    
    # energie per field
    NumberOfEnergies <- rep(0, NumberOfFields)
    for (i in 1:NumberOfFields) {
      NumberOfEnergies[i] <- length(unique(beams$energy[beams$gantryAngle==GantryAngle[i] & beams$patientAngle==PatientSupportAngle[i]]))
    }
    
    
    # write Fluka file
    fluka <- file(paste(file.name, '.fluka', sep=''), "w")
    
    # header
    cat(paste('patient_id ', patient_id, '\n', sep=''), file=fluka, sep='')
    
    # machines (fields?)
    for(F in 1:NumberOfFields) {
      
      beams.machine <- subset(beams, gantryAngle==GantryAngle[F] &
                                patientAngle==PatientSupportAngle[F])
      energies <- unique(beams.machine$energy)
      NumberOfEnergies <- length(energies)
      
      cat(paste('machine# ', F-1, '\n', sep=''), file=fluka, sep='')
      cat(paste('projectile ', ion, '\n', sep=''), file=fluka, sep='')
      cat(paste('charge ', charge, '\n', sep=''), file=fluka, sep='')
      cat(paste('mass ', mass, '\n', sep=''), file=fluka, sep='')
      cat(paste('bolus ', bolus, '\n', sep=''), file=fluka, sep='')
      cat(paste('ripplefilter ', ripplefilter, '\n', sep=''), file=fluka, sep='')
      cat(paste('#submachines ', NumberOfEnergies, '\n', sep=''), file=fluka, sep='')
      cat(paste('#particles ',
                min(beams.machine$fluence),
                max(beams.machine$fluence),
                sum(beams.machine$fluence), 
                '\n', sep=' '), file=fluka, sep='')
      # submachines (layers/energies)
      for(E in 1:NumberOfEnergies) {
        beams.submachine <- subset(beams.machine, energy==energies[E])
        cat(paste('submachine# ',
                  E-1,
                  energies[E],
                  0, # intensità del fascio (non usato),
                  9.4, # FWHM spot
                  '\n', sep=' '), file=fluka, sep='')
        cat(paste('#particles ',
                  min(beams.submachine$fluence),
                  max(beams.submachine$fluence),
                  sum(beams.submachine$fluence), 
                  '\n', sep=' '), file=fluka, sep='')
        cat(paste('stepsize ', stepsize, '\n', sep=''), file=fluka, sep='')
        # points (spots)
        cat(paste('#points ', nrow(beams.submachine), '\n', sep=''), file=fluka, sep='')
        for(S in 1:nrow(beams.submachine)) {
          cat(paste(beams.submachine$deflX[S],
                    beams.submachine$deflY[S],
                    beams.submachine$fluence[S],
                    '\n', sep=' '), file=fluka, sep='')
        }
      }
    }
    
    close(fluka)
  }
}


#' recupera isocentro
#' 
#' @family Beams
#' @export
get.isocenter <- function(plan)
{
  beams <- suppressMessages(get.beams(plan))
  isocenter <- unique(beams[c('x_iso', 'y_iso', 'z_iso')])
  return(isocenter)
}


#' Add field ID column to beams
#' 
#' A unique field ID is associated to a unique combination of gantry and patient angle.
#' 
#' @param beams the beams data.frame
#' 
#' @export
#' @family Beams

add.field <- function(beams)
{
  # check per vedere se esiste già il campo field
  if('field' %in% names(beams)) {
    beams <- beams[, !names(beams)=='field']
  }
  
  beams$field <- interaction(beams[c('gantryAngle', 'patientAngle')])
  levels(beams$field) <- 1:length(levels(beams$field))
  
  return(beams)
}

