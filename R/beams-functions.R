# BEAMS ------------------------------------------------------------------------

#' Create a field object (dataframe)
#'
#' The field object is a dataframe. Different fields are added by apending rows to the dataframe. Each of the following parameter is a dataframe column.
#' @param N number of fields.
#' @param targetVOIIndex the index of the target VOI.
#' @param iecGantryAngle gantry angle.
#' @param iecPatientSupportAngle patient support angle.
#' @param interSpotSpacing spacing among spots along the three axis of the beam port (mm).
#' @param spotsExtensionOutsideTarget thickness of the extension of the spot coverage outside the target VOI.
#' @param targetIsocenter coordinates of the isocenter of the target VOI.
#'
#' @family Beams
#' @export
create.field <- function(N=1,
                         beamLine='protonSimple',
                         targetVOI='PTV',
                         targetVOIIndex=NA,
                         iecGantryAngle=0,
                         iecPatientSupportAngle=0,
                         interSpotSpacing=c(1.2,1.2,1.2),
                         spotsExtensionOutsideTarget=0,
                         targetIsocenter=c(NA,NA,NA)
                         )
{
  field <- data.frame(field=1:N,
                      beamLine=beamLine,
                      targetVOI=targetVOI,
                      targetVOIIndex=targetVOIIndex,
                      iecGantryAngle=iecGantryAngle,
                      iecPatientSupportAngle=iecPatientSupportAngle,
                      interSpotSpacing.x=interSpotSpacing[1],
                      interSpotSpacing.y=interSpotSpacing[2],
                      interSpotSpacing.z=interSpotSpacing[3],
                      spotsExtensionOutsideTarget=spotsExtensionOutsideTarget,
                      targetIsocenter.x=targetIsocenter[1],
                      targetIsocenter.y=targetIsocenter[2],
                      targetIsocenter.z=targetIsocenter[3], stringsAsFactors=FALSE)
  #class(field) <- 'field.plankit'
  field$field <- 1:nrow(field)
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
  beams <- read.table(beams.file, skip=1, stringsAsFactors = FALSE)

  # mette i nomi delle colonne
  if(ncol(beams)==14) {
    names(beams) <- c('particle', 'beamLine',  'x_iso', 'y_iso', 'z_iso', 'gantryAngle', 'patientAngle', 'fluence', 'energy', 'deflX', 'deflY', 'x_s', 'y_s', 'z_s')
    message('Spot positions present...')
  } else if(ncol(beams)==11) {
    names(beams) <- c('particle', 'beamLine', 'x_iso', 'y_iso', 'z_iso', 'gantryAngle', 'patientAngle', 'fluence', 'energy', 'deflX', 'deflY')
    message('Spot positions NOT present...')
  } else if(ncol(beams)==12) {
    names(beams) <- c('x_iso', 'y_iso', 'z_iso', 'gantryAngle', 'patientAngle', 'fluence', 'energy', 'deflX', 'deflY', 'x_s', 'y_s', 'z_s')
    message('Spot positions present (old format)...')
  } else if(ncol(beams)==9) {
    names(beams) <- c('x_iso', 'y_iso', 'z_iso', 'gantryAngle', 'patientAngle', 'fluence', 'energy', 'deflX', 'deflY')
    message('Spot positions NOT present (old format)...')
  }
  #class(beams) <- 'beams'
  return(beams)
}

## Read beams from file (Fluka)
## DA RISCRIVERE PER LA NUOVA Struttura dei beams!!!!
##
## Read RtIonPlan file (Fluka format) and return a beams object.
## @param fluka.file The file name.
## @return a beams object.
## @family Beams
## @export
# read.beams.fluka <- function(fluka.file)
# {
#   message('reading beams file...')
#   beams.fluka <- readLines(fluka.file)
#   index <- 0
#   Nbeams <- 0
# 
#   # ciclo preliminare per numero totale di fasci
#   for(i in 1:length(beams.fluka)) {
#     if(length(grep('#points', beams.fluka[i], fixed=TRUE))>0) {
#       st <- strsplit(beams.fluka[i], ' +')[[1]]
#       Nbeams <- Nbeams + as.numeric(st[2])
#     }
#   }
#   message('number of beams: ', Nbeams)
# 
#   beams <- create.beam(nbeams=Nbeams, with.spots=FALSE)
# 
#   i <- 0
#   while(i < length(beams.fluka)) {
#     i <- i+1
#     if(length(grep('submachine#', beams.fluka[i], fixed=TRUE))>0) {
#       st <- strsplit(beams.fluka[i], ' +')[[1]]
#       energy <- as.numeric(st[3])
#       message(i, ', energy: ', energy)
#     } else if (length(grep('#points', beams.fluka[i], fixed=TRUE))>0) {
#       while(
#         length(grep('submachine#', beams.fluka[i], fixed=TRUE))==0 & i < length(beams.fluka)
#         ) {
#         i <- i+1
#         if(length(grep('submachine#', beams.fluka[i], fixed=TRUE))==0) {
#           st <- strsplit(beams.fluka[i], ' +')[[1]]
#           deflX <- as.numeric(st[1])
#           deflY <- as.numeric(st[2])
#           fluence <- as.numeric(st[3])
#           index <- index+1
#           beams$energy[index] <- energy
#           beams$deflX[index] <- deflX
#           beams$deflY[index] <- deflY
#           beams$fluence[index] <- fluence
#           if(is.na(deflX)) {message('NA: ', i)}
#         }
#         #print(beams[index,])
#       }
#       i <- i-1
#     }
#   }
# 
#   #class(beams) <- beams
#   return(beams)
# }


#' Get beams from plan
#'
#' Get the beams object associated to the plan. There are two type of beams: input and output. Inputs beams indicates the beams used for a forward planning, or for the initial condition of the optimization procedure in the inverse planning. Output beams indicates the beams evaluated in the inverse planning.
#' @param plan The plan object.
#' @param input If TRUE the input beams are retrieved.
#' @family Beams
#' @export
get.beams <- function(plan, ...) UseMethod('get.beams')

#' @family Beams
#' @export
get.beams.plankit.plan <- function(plan, input=FALSE)
{
  if(!is.null(plan$beams)) { # rimane qui un problema di identificazione di cosa sono i beams...
    return(plan$beams)
  } else {
    if(input){beams.file <- get.filepath('inputBeamsFile', plan=plan)}
    else {beams.file <- get.filepath('outputBeamsFile', plan=plan)}
    if(is.null(beams.file)) {
      stop('beams object not present')
    } else {
      return(read.beams(beams.file))
    }
  }
}

#' @family Beams
#' @export
get.beams.gate.plan <- function(plan)
{
  if(!is.null(plan$beams)) {
    return(plan$beams)
  } else {
    beams.file <- paste(plan$name, plan$beamsFile, sep='/')
    return(read.beams(beams.file))
  }

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


#' Create an empty beam(s) object
#'
#' The beam(s) object is a dataframe with columns:
#' \itemize{
#'  \item{"x_iso, y_iso, z_iso"}{Stuff}
#'  \item{"parameter 2"}{Stuff}
#' }
#' @param nbeams The number of the beams
#' @param particle particle name
#' @param beamLine the beamLine name (i.e. the corresponding LUT to use).
#' @param x_iso,y_iso,z_iso Coordinates of the isocenter (mm).
#' @param gantryAngle,patientAngle Gantry and patient support angles (degree).
#' @param fluence The number of particles for the beam(s).
#' @param energy The kinetic energy of the particles (MeV/u).
#' @param deflX,deflY Deflections on the plane normal to the beam direction and passing through the isocenter (mm).
#' @param x_s,y_s,z_s Spot coordinates, i.e. the position of the Bragg peak (mm).
#' @family Beams
#' @export
create.beam <- function(particle='1H', beamLine='simpleProtons' ,x_iso=0, y_iso=0, z_iso=0, gantryAngle=0, patientAngle=0, fluence=1, energy=0, deflX=0, deflY=0, x_s=0, y_s=0, z_s=0)
{
  beam <- data.frame(particle=particle, beamLine=beamLine, x_iso=x_iso,  y_iso=y_iso, z_iso=z_iso, gantryAngle=gantryAngle, patientAngle=patientAngle, fluence=fluence, energy=energy, deflX=deflX, deflY=deflY, x_s=x_s, y_s=y_s, z_s=z_s)
  return(beam)
}


#' Write beams (RTIonPlan)
#'
#' Write the beams object on file. It can use different formats (default PlanKIT).
#'
#' @param beams The beams object
#' @param file.name The file name.
#' @param format The format. It can be 'plankit', 'gate' or 'fluka'.
#' @param ion. The primary ion. (used by fluka).
#' @family Beams
#' @export
#'
write.beams <- function(beams, file.name, format='puredek', ion='1H', add.extension=TRUE)
{

  # PlanKIT
  if(format=='plankit' | format=='puredek') {
    if(add.extension) {file.name <- paste(file.name, '.beams', sep='')}
    N <- nrow(beams)
    con <- file(file.name, "w")
    writeLines(paste(N), con=con)
    close(con)

    # elimina colonna fields
    index.field <- which(names(beams)=='field')
    if(length(index.field>0)) {beams <- beams[,-index.field]}

    # riordina le colonne per essere consistente con pure-dek
    # spot position presenti

    if(all(c("particle", "beamLine", "x_iso", "y_iso", "z_iso", "gantryAngle", "patientAngle", "fluence", "energy", "deflX", "deflY", "x_s", "y_s", "z_s") %in% colnames(beams))) {
      beams.ord <- beams[c("particle", "beamLine", "x_iso", "y_iso", "z_iso", "gantryAngle", "patientAngle", "fluence", "energy", "deflX", "deflY", "x_s", "y_s", "z_s")]
    } else if(all(c("particle", "beamLine", "x_iso", "y_iso", "z_iso", "gantryAngle", "patientAngle", "fluence", "energy", "deflX", "deflY") %in% colnames(beams))) {
      beams.ord <- beams[c("particle", "beamLine", "x_iso", "y_iso", "z_iso", "gantryAngle", "patientAngle", "fluence", "energy", "deflX", "deflY")]
    } else {
      stop('Error, missing mandatory data in beams.')
    }

    write.table(beams.ord, file=file.name, append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
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
    if(add.extension) {file.name <- paste(file.name, '.gate', sep='')}
    gate <- file(file.name, "w")

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
    if(add.extension) {file.name <- paste(file.name, '.fluka', sep='')}
    fluka <- file(file.name, "w")

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
                  E, # E-1.
                  energies[E],
                  1, #0, # intensità del fascio (non usato),
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
#' An unique field ID is associated to a unique combination of gantry and patient angle.
#'
#' @param beams the beams data.frame
#' @param numeric use a numeric ID for the field
#'
#' @export
#' @family Beams
add.field <- function(beams, numeric=TRUE)
{
  # check per vedere se esiste già il campo field
  if('field' %in% names(beams)) {
    beams <- beams[, !names(beams)=='field']
  }

  beams$field <- interaction(beams[c('beamLine', 'gantryAngle', 'patientAngle')])
  if(numeric) {
    levels(beams$field) <- 1:length(levels(beams$field))
    beams$field <- as.factor(as.numeric(beams$field))
  }

  return(beams)
}

#' Get rays from beams
#'
#' Get rays (normalized vector directions + absolute coordinate) in the reference frame of the CT.
#'
#' @param beams the beams data.frame
#' @param s distance of the source from the isocenter (it uses only one virtual source).
#' @param unique.rays Return only unique rays
#' @return data frame of rays. Each ray has 6 component: a point (X,Y,Z) on the plane passing through the isocenter, and a normalized vector for the direction (xn, yn, zn).
#' @export
#' @family Beams
get.rays <- function(beams,
                     s=4000,
                     unique.rays=FALSE)
{

  X0 <- X1 <- X2 <- c(1,0,0)
  Y0 <- Y1 <- Y2 <- c(0,0,1) # = Z nel sistema di riferimento della CT
  s0 <- s1 <- s2 <- c(0,-s,0)

  # per arrivare nel sistema di riferimento della CT occorre ruotare tutto di 90
  # attorno all'asse X

  # il gantry corrisponde a ruotare attorno all'asse Z del sistema di riferimento della CT
  # il patient angel corrisponde a ruotare attorno all'asse X (?)

  Nb <- nrow(beams)

  rays <- data.frame(X=rep(0, Nb), Y=0, Z=0, xn=0, yn=0, zn=0)

  for(i in 1:Nb) {
    gantryAngle <- beams$gantryAngle[i]
    patientAngle <- beams$patientAngle[i]

    # gantry (attorno a z)
    X1[1] <- (X0[1]*cos(gantryAngle) - X0[2]*sin(gantryAngle))
    X1[2] <- (X0[1]*sin(gantryAngle) + X0[2]*cos(gantryAngle))
    Y1[1] <- (Y0[1]*cos(gantryAngle) - Y0[2]*sin(gantryAngle))
    Y1[2] <- (Y0[1]*sin(gantryAngle) + Y0[2]*cos(gantryAngle))
    s1[1] <- (s0[1]*cos(gantryAngle) - s0[2]*sin(gantryAngle))
    s1[2] <- (s0[1]*sin(gantryAngle) + s0[2]*cos(gantryAngle))

    #patient (attorno a x)
    X2[1] <- X1[1]
    X2[3] <- (X1[3]*cos(patientAngle) - X1[2]*sin(patientAngle))
    X2[2] <- (X1[3]*sin(patientAngle) + X1[2]*cos(patientAngle))
    Y2[1] <- Y1[1]
    Y2[3] <- (Y1[3]*cos(patientAngle) - Y1[2]*sin(patientAngle))
    Y2[2] <- (Y1[3]*sin(patientAngle) + Y1[2]*cos(patientAngle))
    s2[1] <- s1[1]
    s2[3] <- (s1[3]*cos(patientAngle) - s1[2]*sin(patientAngle))
    s2[2] <- (s1[3]*sin(patientAngle) + s1[2]*cos(patientAngle))

    R <- X2*beams$deflX[i]+Y2*beams$deflY[i] + c(beams$x_iso[i], beams$y_iso[i], beams$z_iso[i])
    rn <- R-s2
    rn <- rn / sqrt(sum(rn^2))

    rays$X[i] <- R[1]
    rays$Y[i] <- R[2]
    rays$Z[i] <- R[3]
    rays$xn[i] <- rn[1]
    rays$yn[i] <- rn[2]
    rays$zn[i] <- rn[3]

    #X <- R*X0 * beams$deflX[i]
    #Y <- R*Y0 * beams$deflY[i]
  }

  if(unique.rays) {
    rays <- rays[!duplicated(rays),]
  }

  return(rays)
}


#' Get entrance point in CT
#'
#' Get the coordinates of the entrance point in the CT for the specified rays.
#'
#' @param rays the rays data.frame
#' @return data frame of rays in wich the coordinates of the entrance point are added (Xin, Yin, Zin).
#' @export
#' @family Beams
get.ct.entrance.point <- function(rays, ct)
{

  # coordinate piani esterni
  xi <- range(create.intervals(ct$x))
  yi <- range(create.intervals(ct$y))
  zi <- range(create.intervals(ct$z))

  for(ii in 1:nrow(rays)) {
    ray <- rays[ii,]
    p <- ray.tracing(ray, xi, yi, zi)

    # trova ingresso
    rays$Xin[ii] <- p$x[1]
    rays$Yin[ii] <- p$y[1]
    rays$Zin[ii] <- p$z[1]
  }

  return(rays)
}
