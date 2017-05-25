# mette tutte le funzioni di visualizzazione

# LUT DISPLAY FUNCTIONS --------------------------------------------------------

#' 2D display of the beam-lut
#'
#' The \code{lut} object use a structure similar to the "values" object. For example the different quantities are stored in a 5D array, \code{lut$values[i,j,k,l]}. The indices of the array correspond to: variable, energy, x, y, normalized z.
#'
#' @param lut the Look Up Table object
#' @param variable the specific variable to be displayed (not needed if values contains only one variable, or if the variable is the first of the list)
#' @param E the energy of the beam (if it is not specified, the mean energy of the available energies is used). The energy displayed is approximated to the nearest available in the lut.
#' @param bw display using gray color levels (optional, boolean)
#' @param r.cut cut-off radius (optional)
#' @param cint display contour levels (optional, boolean)
#' @param z.absolute display absolute z (optional, boolean)
#' @export
display.lut <- function(lut, variable=NULL, E=NULL, bw=FALSE, r.cut=NULL, cont=FALSE, z.absolute=FALSE)
{
  #suppressMessages(library(fields))

  #cat(dim(lut$values), '\n')

  # variable
  v <- NULL
  if(!is.null(variable)) {v <- which(lut$variables==variable)}
  if(length(v)==0) {v <- 1}
  variable <- lut$variable[v]
  message('displaying ', variable, ' (', v, ')...')
  lut$values <- lut$values[v,,,,,drop=FALSE]
  #cat(dim(lut$values), '\n')

  # energy
  if(is.null(E)) {E <- mean(lut$E)}
  message('plotting E=',E)
  NE <- which(abs(lut$E-E) == min(abs(lut$E-E)))[1]
  lut$values <- lut$values[1,NE,,,]
  #cat(dim(lut$values), '\n')

  # plane
  if(lut$Ny==1) {
    message('radial lut')
  } else {
    message('plotting middle plane (y=0)...')
    Ny <- which(abs(lut$y) == min(abs(lut$y)))[1]
    lut$values <- lut$values[,Ny,]
  }

  # cut-off
  if(!is.null(r.cut)) {
    message('using r.cut=',r.cut)
    index.x <- abs(lut$x) < r.cut
    lut$values <- lut$values[index.x,]
    lut$x <- lut$x[index.x]
  }

  # Z
  if(z.absolute) {
    Z <- lut$z * lut$z.BP(E)
  } else {
    Z <- lut$z
  }

  # colori
  Nc <- 255
  if(bw) {
    cols <- gray( ((1:Nc)-1)/Nc )
  } else {
    cols <- rainbow(Nc, start=0, end=0.6)[Nc:1]
  }

  # plot immagine
  main <- paste('E = ', E, ' MeV/u\n', variable, sep='')
  subt <- ''
  if(lut$Ny==1) {my.xlab <- 'r [mm]'} else {my.xlab <- 'x [mm]'}
  if(z.absolute) {my.ylab <- 'z [mm]'} else {my.ylab <- 'z/zBp'}
  fields::image.plot(lut$x, Z, lut$values, col=cols,
             main=main, sub=subt, xlab=my.xlab, ylab=my.ylab)

  # plot contorni
  if(cont==TRUE) {
    contour(lut$x, Z, lut$values, add=TRUE)
  }
}


# VALUES DISPLAY FUNCTIONS -----------------------------------------------------


#' Display a 2D slice
#'
#' Display an arbitrary 2D slice (axial, coronal, sagittal) of the values stored in a \code{values} object.
#'
#' If slice index and slice coordinate are not specified, but the plan is specified, the displayed slice is axial and at the isocenter of the plan.
#' If slice index, slice coordinate and plan are not specified, the displayed slice is axial and at the middle of the z range.
#'
#' If the plan is specified, it is possible to not specify the values object. In this case the one associated to the plan is used.
#'
#' @param values the values object (optional if plan is defined)
#' @param variable the specific variable to be displayed (not needed if values contains only one variable, or if the variable is the first of the list)
#' @param Nx,Ny,Nz slice index to be displayed. Only one of them should be defined, corresponding to an axial, coronal or sagittal slice (optional)
#' @param x,y,z slice coordinate to be displayed. Only one of them should be defined, corresponding to an axial, coronal or sagittal slice (optional)
#' @param plan the plan object (optional if values is defined)
#' @param cont display isolevel contours of the values (optional, boolean)
#' @param isoc display the isocenter on the slice (optional, boolean)
#' @param bw display using gray color levels. It colud be used to display CT values (optional, boolean)
#' @param xlim,ylim,zlim two element vectors containing the coordinate ranges to be displayed (optional)
#' @param vlim two element vector containing the value range to be displayed. Values outside the range are displayed as 'blank', i.e. no color tile (optional)
#' @param colors vector of color values (colormap)
#' @param file.name name of the file on which to save the figure. By default no file is written (optional)
#' @param width,height sizes of the figure to be saved (inches) (optional)
#'
#' @return used slice (index and coordinates) data.frame. It is used to make further combination with other plots and procedures (e.g. \code{display.slice.all})
#'
#' @family display slices
#' @export
#' @importFrom fields set.panel image.plot

display.slice <- function(values=NULL,
                          variable=NULL,
                          Nx=NA, Ny=NA, Nz=NA,
                          x=NA, y=NA, z=NA,
                          plan=NULL,
                          bw=FALSE,
                          cont=FALSE,
                          isoc=FALSE,
                          xlim=NULL,
                          ylim=NULL,
                          zlim=NULL,
                          vlim=NULL,
                          colors=NULL,
                          invert.y.axis=FALSE,
                          file.name=NULL,
                          width=7, height=7,
                          levels=NULL)
{
  #suppressMessages(library(fields))

  # recupera values da plan se values=NULLs
  if(is.null(values)) {values <- get.values(plan)}

  # variable
  v <- NULL
  if(!is.null(variable)) {v <- which(values$variables==variable)}
  if(length(v)==0) {v <- 1}
  variable <- values$variable[v]
  message('displaying ', variable, ' ...')
  if(values$Nv>1) {values$values <- values$values[v,,,]}

  # recupera isocentro
  isocenter <- NULL
  Nx.iso <- Ny.iso <- Nz.iso <- NA
  ISOCENTER <- FALSE
  if(!is.null(plan)) {
    isocenter <- get.isocenter(plan)
    x_iso <- isocenter$x_iso
    y_iso <- isocenter$y_iso
    z_iso <- isocenter$z_iso
    Nx.iso <- identify.slice(x_iso, values$x)
    Ny.iso <- identify.slice(y_iso, values$y)
    Nz.iso <- identify.slice(z_iso, values$z)

    message('isocenter: (', x_iso, ",", y_iso, ",", z_iso, ')')
    message('isocenter voxel: (', Nx.iso, ",", Ny.iso, ",", Nz.iso, ')')
  }

  # estremi
  if(is.null(xlim)) {xlim <- range(values$x, na.rm=TRUE)}
  if(is.null(ylim)) {ylim <- range(values$y, na.rm=TRUE)}
  if(is.null(zlim)) {zlim <- range(values$z, na.rm=TRUE)}
  if(is.null(vlim)) {vlim <- range(values$values, na.rm=TRUE)}
  if(invert.y.axis) {ylim <- rev(ylim)}

  # slice...

  # coordinate specificate
  if(!is.na(x)) {Nx <- identify.slice(x, values$x)}
  if(!is.na(y)) {Ny <- identify.slice(y, values$y)}
  if(!is.na(z)) {Nz <- identify.slice(z, values$z)}

  # coordinate non specificate, indici specificati
  if(!is.na(Nx) & is.na(x)) {x <- values$x[Nx]}
  if(!is.na(Ny) & is.na(y)) {y <- values$y[Ny]}
  if(!is.na(Nz) & is.na(z)) {z <- values$z[Nz]}

  # niente di specificato
  if(is.na(Nx) & is.na(Ny) & is.na(Nz)
     & is.na(x) & is.na(y) & is.na(z)) {
    if(!is.null(plan)) {
      # isocentro
      z <- z_iso
      Nz <- Nz.iso
      ISOCENTER <- TRUE
    } else {
      Nz <- round(values$Nz/2)
      z <- values$z[Nz]
    }
  }

  if(!is.null(isocenter)) {
    if(!is.na(Nx)) {if(Nx==Nx.iso) ISOCENTER <- TRUE}
    if(!is.na(Ny)) {if(Ny==Ny.iso) ISOCENTER <- TRUE}
    if(!is.na(Nz)) {if(Nz==Nz.iso) ISOCENTER <- TRUE}
  }

  message('plotting voxels: (', Nx, ",", Ny, ",", Nz, ')')
  message('plotting coord.: (', x, ",", y, ",", z, ')')

  # text
  if(!is.null(plan)) {main <- paste(plan$name, '\n', variable, sep='')} else {main <- variable}
  if(!is.na(z)) {subt <- paste('slice at z =', z, 'mm')}
  else if(!is.na(y)) {subt <- paste('slice at y =', y, 'mm')}
  else if(!is.na(x)) {subt <- paste('slice at x =', x, 'mm')}


  # colori
  if(!is.null(colors)) {
    cols <- colors
  } else {
    Nc <- 255
    if(bw) {
      cols <- gray( ((1:Nc)-1)/Nc )
    } else {
      cols <- rainbow(Nc, start=0, end=0.6)[Nc:1]
    }
  }

  # INIZIA FIGURA
  if(!is.null(file.name)) {
    png(filename=file.name, width=width, height=height, res=300, units='in')
  }

  # panel
  set.panel()
  par(oma=c(0,0,0,4)) # margin of 4 spaces width at right hand side
  set.panel(1,1) # 1X1 matrix of plots

  # plot immagine
  if(!is.na(Nz)) {
    image.plot(values$x, values$y, as.matrix(values$values[,,Nz]), col=cols,
               xlim=xlim, ylim=ylim, zlim=vlim,
               main=main, sub=subt, xlab='x [mm]', ylab='y [mm]')
  } else if(!is.na(Ny)) {
    image.plot(values$x, values$z, as.matrix(values$values[,Ny,]), col=cols,
               xlim=xlim, ylim=zlim, zlim=vlim,
               main=main, sub=subt, xlab='x [mm]', ylab='z [mm]')
  } else if(!is.na(Nx)) {
    image.plot(values$y, values$z, as.matrix(values$values[Nx,,]), col=cols,
               xlim=ylim, ylim=zlim, zlim=vlim,
               main=main, sub=subt, xlab='y [mm]', ylab='z [mm]')
  }

  # plot contorni
  if(is.null(levels)) {levels <- pretty(vlim, 10)}
  if(cont==TRUE) {
    if(!is.na(Nz)) {contour(values$x, values$y, values$values[,,Nz], add=TRUE, levels=levels)}
    else if(!is.na(Ny)) {contour(values$x, values$z, values$values[,Ny,], add=TRUE, levels=levels)}
    else if(!is.na(Nx)) {contour(values$y, values$z, values$values[Nx,,], add=TRUE, levels=levels)}
  }


  # plot isocentro
  if(isoc) {
    points(isocenter$x_iso, isocenter$y_iso)
    if(ISOCENTER) {
      abline(v=isocenter$x_iso, lty=2)
      abline(h=isocenter$y_iso, lty=2)
    }
  }

  par(oma=c( 0,0,0,1))# reset margin to be much smaller.
  #image.plot(legend.only=TRUE, zlim=vlim, col=cols)

  set.panel()

  # FINISCE FIGURA
  if(!is.null(file.name)) {dev.off(); message('figure saved in ', file.name)}

  # ritorna informazioni sulla slice
  return(data.frame(Nx=Nx, Ny=Ny, Nz=Nz, x=x, y=y, z=z))
}


#' Display CT slices
#'
#' Alias for display.slices().
#' @family display slices
#' @export
#' @importFrom fields set.panel image.plot
display.slice.ct <- function(ct, contours=NULL,
                             x=NA, y=NA, z=NA,
                             Nx=NA, Ny=NA, Nz=NA,
                             plan=NULL,
                             xlim=NULL,
                             ylim=NULL,
                             zlim=NULL,
                             vlim=NULL,
                             HU.window=c(-1000,3000),
                             invert.y.axis=FALSE,
                             file.name=NULL,
                             width=7, height=7, dpi=300,
                             use.contours.colors=TRUE,
                             contours.legend=FALSE,
                             cex.contours.legend=0.8,
                             contour.legend.position='topleft')
{

  # estremi
  if(is.null(xlim)) {xlim <- range(ct$x)}
  if(is.null(ylim)) {ylim <- range(ct$y)}
  if(is.null(zlim)) {zlim <- range(ct$z)}
  if(is.null(vlim)) {vlim <- range(ct$values, na.rm=TRUE)}
  if(invert.y.axis) {ylim <- rev(ylim)}
  #vlim <- range(ct$values, na.rm=TRUE)

  
  # colori
  my.colors <- colormap.ct(HU.range=vlim, HU.window=HU.window)
  if(use.contours.colors | contours.legend) {
    if( !('display.color' %in% colnames(contours)) ) {
      contours <- add.colours.contours(contours)
    }
    use.contour.colors <- TRUE # se si vuole visualizzare la legenda allora occorre usare i colori...
  }

  # INIZIA FIGURA
  if(!is.null(file.name)) {
    png(filename=file.name, width=width, height=height, res=dpi, units='in')
  }

  # immagine ct
  slice <- display.slice(ct, variable='HounsfieldNumber',
                         Nx=Nx, Ny=Ny, Nz=Nz,
                         x=x, y=y, z=z,
                         plan=plan, bw=TRUE,
                         xlim=xlim,
                         ylim=ylim,
                         zlim=zlim,
                         vlim=vlim,
                         colors=my.colors)
  Nx <- slice$Nx
  Ny <- slice$Ny
  Nz <- slice$Nz

  # panel
  set.panel()
  par(oma=c(0,0,0,4)) # margin of 4 spaces width at right hand side
  set.panel(1,1) # 1X1 matrix of plots

  # contorni (solo per slice assiali)
  if(!is.null(contours) & !is.na(Nz)) {
    message('using contours...')
    roi.s <- subset(contours, slice==(Nz-1)) # nota: il numero della slice dei contorni inizia da zero, mentre l'indice-slice della CT da 1.
    roi.names <- unique(roi.s$contour)
    N.roi <- length(roi.names)
    if(N.roi>0) {
    message('adding contours for: ', paste(roi.names, collapse=' '))
    my.lwd <- rep(2, N.roi)
    for(i in 1:N.roi) {
      rs <- subset(roi.s, contour==roi.names[i])
      type <- unique(rs$type)
      if(type=='PTV') {my.lwd[i] <- 4}
      pol <- unique(rs$polygon)
      if(use.contour.colors) { # imposta colori predefiniti
        c.col <- rs$display.color[1]
      } else {c.col <- 'green'}
      #print(pol)
      for(j in 1:length(pol)) {
        rsp <- subset(rs, polygon==pol[j])
        rsp <- rbind(rsp, rsp[1,])
        lines(rsp$x, rsp$y, col=c.col, lwd=my.lwd[i])
      }
    }
    }
  }

  # legenda contorni
  if(contours.legend) {
    display.contours.legend(contours = contours, position = contour.legend.position, cex = cex.contours.legend)
  }
  
  par(oma=c(0,0,0,1))# reset margin to be much smaller.
  #image.plot(legend.only=TRUE, zlim=vlim, col=col.val)

  set.panel()

  
  
  # FINISCE FIGURA
  if(!is.null(file.name)) {dev.off(); message('figure saved in ', file.name)}

  return(slice)
}


#' Display a 2D slice superimposed to the CT
#'
#' Display an arbitrary 2D slice (axial, coronal, sagittal) of the values stored in a \code{values} object. The slice is superimposed to the associated CT slice. Optionally ROI contours can be also displayed.
#'
#' If slice index and slice coordinate are not specified, but the plan is specified, the displayed slice is axial and at the isocenter of the plan.
#' If slice index, slice coordinate and plan are not specified, the displayed slice is axial and at the middle of the z range.
#'
#' If the plan is specified, it is possible to not specify some or all of the following objects: values, CT and contours. In this case the missing ones are replaced with those associated to the plan.
#'
#' @param ct the ct object (optional, if plan is specified)
#' @param contours the contours object (optional)
#' @param values the values object (optional, if plan is specified)
#' @param variable the specific variable to be displayed (not needed if values contains only one variable, or if the variable is the first of the list)
#' @param Nx,Ny,Nz slice index to be displayed. Only one of them should be defined, corresponding to an axial, coronal or sagittal slice (optional)
#' @param x,y,z slice coordinate to be displayed. Only one of them should be defined, corresponding to an axial, coronal or sagittal slice (optional)
#' @param plan the plan object (optional, if \code{values} and \code{ct} are defined)
#' @param cont display isolevel contours of the values (optional, boolean)
#' @param isoc display the isocenter on the slice (optional, boolean)
#' @param bw display using gray color levels. It colud be used to display CT values (optional, boolean)
#' @param xlim,ylim,zlim two element vectors containing the coordinate ranges to be displayed (optional)
#' @param vlim two element vector containing the value range to be displayed. Values outside the range are displayed as 'blank', i.e. no color tile (optional)
#' @param file.name name of the file on which to save the figure. By default no file is written (optional)
#' @param width,height sizes of the figure to be saved (inches) (optional)
#' @param HU.window HU window to visualize for the CT.
#' @param alpha.lower,alpha.upperThe opacity (alpha) is evaluated as a linear ramp (\code{alpha.lower} to \code{alpha.upper}) from \code{alpha.lower} to \code{alpha.upper} of the values range.
#' @param invert.y.axis Invert the y axis.
#' @param contour.color The color used for the contour of the VOI.
#' @param dpi dpi of the saved image.
#' @param use.contour.colors Use the colors defined in the contours object for the contours.
#' @param display.contours.legend Display a legend for the contours.
#' @param cex.contours.legend Font dimension for the contour legend.
#' @param contour.legend.position Position of the contour legend.
#' @param title Explicitly define the text to be displayed in the plot title.
#' @param col.invert Invert the color map (red for lower values and blue for high values)
#'
#' @family display slices
#' @export
#' @importFrom fields set.panel image.plot

display.slice.all <- function(ct=NULL,
                              contours=NULL,
                              values=NULL,
                              variable=NULL,
                              plan=NULL,
                              x=NA, y=NA, z=NA,
                              cont=FALSE,
                              isoc=FALSE,
                              xlim=NULL,
                              ylim=NULL,
                              zlim=NULL,
                              vlim=NULL,
                              file.name=NULL,
                              width=7, height=7,
                              HU.window=c(-1000,3000),
                              alpha.lower=0,
                              alpha.upper=1,
                              invert.y.axis=FALSE,
                              contour.color='green',
                              dpi=300,
                              levels=NULL,
                              use.contour.colors=TRUE,
                              contours.legend=FALSE,
                              cex.contours.legend=0.8,
                              contour.legend.position='topleft',
                              title=NULL,
                              col.invert=FALSE)
{
  #suppressMessages(library(fields))

  # recupera oggetti da plan
  if(is.null(values)) {message(class(plan)); values <- get.values(plan)}
  if(is.null(ct)) {ct <- get.ct(plan)}
  if(is.null(contours) & !is.null(plan)) {contours <- get.contours(plan)}
  #if(!use.contours) {contours <- NULL}

  Nx.ct <- Ny.ct <- Nz.ct <- Nx.values <- Ny.values <- Nz.values <- NA

  # recupera isocentro
  isocenter <- NULL
  x_iso <- y_iso <- x_iso <- NA
  ISOCENTER <- FALSE
  if(!is.null(plan)) {
    isocenter <- get.isocenter(plan)
    x_iso <- isocenter$x_iso
    y_iso <- isocenter$y_iso
    z_iso <- isocenter$z_iso
    message('isocenter: (', x_iso, ",", y_iso, ",", z_iso, ')')
  }

  # coordinate specificate
  if(!is.na(x)) {
    Nx.values <- identify.slice(x, values$x)
    Nx.ct <- identify.slice(x, ct$x)
  }
  if(!is.na(y)) {
    Ny.values <- identify.slice(y, values$y)
    Ny.ct <- identify.slice(y, ct$y)
  }
  if(!is.na(z)) {
    Nz.values <- identify.slice(z, values$z)
    Nz.ct <- identify.slice(z, ct$z)
  }

  # niente di specificato
  if(is.na(x) & is.na(y) & is.na(z)) {
    if(!is.null(plan)) {
      # isocentro
      z <- z_iso
      ISOCENTER <- TRUE
    } else {
      z <- mean(values$z)
    }
    Nz.values <- identify.slice(z, values$z)
    Nz.ct  <- identify.slice(z, ct$z)
  }

  if(!is.null(isocenter)) {
    #if(!is.na(x_iso)) {if(x==x_iso) ISOCENTER <- TRUE}
    #if(!is.na(y_iso)) {if(y==y_iso) ISOCENTER <- TRUE}
    #if(!is.na(z_iso)) {if(z==z_iso) ISOCENTER <- TRUE}
  }

  message('plotting at coord.: (', x, ",", y, ",", z, ')')
  message('ct voxels: (', Nx.ct, ",", Ny.ct, ",", Nz.ct, ')')
  message('values voxels: (', Nx.ct, ",", Ny.ct, ",", Nz.ct, ')')

  # variable
  v <- which(values$variables==variable)
  if(length(v)==0) {v <- 1}
  variable <- values$variable[v]
  message('displaying ', variable, '...')
  if(values$Nv>1) {values$values <- values$values[v,,,]}

  # text
  if(!is.null(title)) {main <- title} else {
    if(!is.null(plan)) {main <- paste(plan$name, '\n', variable, sep='')} else {main <- variable}
  }
  if(!is.na(z)) {subt <- paste('slice at z =', z, 'mm')}
  else if(!is.na(y)) {subt <- paste('slice at y =', y, 'mm')}
  else if(!is.na(x)) {subt <- paste('slice at x =', x, 'mm')}

  # colori
  Nc <- 255
  Nc.alpha <- round(Nc*(alpha.upper-alpha.lower))
  interval <- ((1:Nc)-1)/Nc
  interval.alpha <- c( rep(alpha.lower,round(Nc*alpha.lower)), seq(alpha.lower, alpha.upper, length.out=Nc.alpha), rep(alpha.upper,round(Nc-Nc*alpha.upper)) )
  col.val <- hsv( interval[Nc:1]*0.64, alpha=interval.alpha[1:Nc])
  col.ct <- colormap.ct(HU.range=range(ct$values), HU.window=HU.window)
  
  if(col.invert) {col.val <- col.val[length(col.val):1]}
  
  # colori contorni
  if(use.contour.colors | contours.legend) {
    if( !('display.color' %in% colnames(contours)) ) {
      contours <- add.colours.contours(contours)
    }
    use.contour.colors <- TRUE # se si vuole visualizzare la legenda allora occorre usare i colori...
  }


  # estremi
  if(is.null(xlim)) {xlim <- range(values$x)}
  if(is.null(ylim)) {ylim <- range(values$y)}
  if(is.null(zlim)) {zlim <- range(values$z)}
  if(is.null(vlim)) {vlim <- range(values$values, na.rm=TRUE)}

  if(invert.y.axis) {ylim <- rev(ylim)}

  # INIZIA FIGURA
  if(!is.null(file.name)) {
    png(filename=file.name, width=width, height=height, res=dpi, units='in')
  }

  set.panel()
  par(oma=c(0,0,0,4)) # margin of 4 spaces width at right hand side
  set.panel(1,1) # 1X1 matrix of plots

  #asp <- 1/abs( diff(range(values$x)) /  diff(range(values$y)) )
  #print(asp)

  # layer colorbar
  if(!is.na(Nz.values)) {
    image(values$x, values$y, values$values[,,Nz.values], col=col.val,
               xlim=xlim, ylim=ylim, zlim=vlim,
               main=main, sub=subt, xlab='x [mm]', ylab='y [mm]')
  } else if(!is.na(Ny.values)) {
    image(values$x, values$z, values$values[,Ny.values,], col=col.val,
               xlim=xlim, ylim=zlim, zlim=vlim,
               main=main, sub=subt, xlab='x [mm]', ylab='z [mm]')
  } else if(!is.na(Nx.values)) {
    image(values$y, values$z, values$values[Nx.values,,], col=col.val,
               xlim=ylim, ylim=zlim, zlim=vlim,
               main=main, sub=subt, xlab='y [mm]', ylab='z [mm]')
  }

  # layer ct
  if(!is.na(Nz.ct)) {
    image(ct$x, ct$y, as.matrix(ct$values[,,Nz.ct]), col=col.ct, xlim=xlim, ylim=ylim, add=TRUE)
  } else if(!is.na(Ny.ct)) {
    image(ct$x, ct$z, as.matrix(ct$values[,Ny.ct,]), col=col.ct, xlim=xlim, ylim=zlim, add=TRUE)
  } else if(!is.na(Nx.ct)) {
    image(ct$y, ct$z, as.matrix(ct$values[Nx.ct,,]), col=col.ct, xlim=ylim, ylim=zlim, add=TRUE)
  }

  # layer values
  if(!is.na(Nz.values)) {
    image(values$x, values$y, values$values[,,Nz.values], col=col.val, add=TRUE, xlim=xlim, ylim=ylim, zlim=vlim)
  } else if(!is.na(Ny.values)) {
    image(values$x, values$z, values$values[,Ny.values,], col=col.val, add=TRUE, xlim=xlim, ylim=zlim, zlim=vlim)
  } else if(!is.na(Nx.values)) {
    image(values$y, values$z, values$values[Nx.values,,], col=col.val, add=TRUE, xlim=ylim, ylim=zlim, zlim=vlim)
  }


  # layer contorni
  if(is.null(levels)) {levels <- pretty(vlim, n=10)}
  if(cont==TRUE) {
    if(!is.na(Nz.values)) {contour(values$x, values$y, values$values[,,Nz.values], add=TRUE, xlim=xlim, ylim=ylim, zlim=vlim, levels=levels)}
    else if(!is.na(Ny.values)) {contour(values$x, values$z, values$values[,Ny.values,], add=TRUE, xlim=xlim, ylim=zlim, zlim=vlim, levels=levels)}
    else if(!is.na(Nx.values)) {contour(values$y, values$z, values$values[Nx.values,,], add=TRUE, xlim=ylim, ylim=zlim, zlim=vlim, levels=levels)}
  }

  # layer roi (solo slice assiali)
  if(!is.null(contours) & !is.na(Nz.values)) {
    message('using voi contours...')
    roi.s <- subset(contours, slice==(Nz.ct-1)) # nota: il numero della slice dei contorni inizia da zero, mentre l'indice-slice della CT da 1.
    roi.names <- unique(roi.s$contour)
    N.roi <- length(roi.names)
    if(N.roi>0) {
    message('adding contours for: ', paste(roi.names, collapse=' '))
    my.lwd <- rep(2, N.roi)
    for(i in 1:N.roi) {
      rs <- subset(roi.s, contour==roi.names[i])
      rs <- rbind(rs, rs[1,]) # chiudi il contorno...
      type <- unique(rs$type)
      if(type=='PTV') {my.lwd[i] <- 4}
      pol <- unique(rs$polygon)
      if(use.contour.colors) { # imposta colori predefiniti
        c.col <- rs$display.color[1]
      } else {c.col <- contour.color}
      for(j in 1:length(pol)) {
        rsp <- subset(rs, polygon==pol[j])
        lines(rsp$x, rsp$y, col=c.col, lwd=my.lwd[i])
      }
    }
    }
  }

  # layer isocentro
  if(isoc) {
    points(isocenter$x_iso, isocenter$y_iso)
    if(ISOCENTER) {
      abline(v=isocenter$x_iso, lty=2)
      abline(h=isocenter$y_iso, lty=2)
    }
  }

  par(oma=c( 0,0,0,1))# reset margin to be much smaller.
  fields::image.plot(legend.only=TRUE, zlim=vlim, col=col.val)
  
  # legenda contorni
  if(contours.legend) {
    display.contours.legend(contours = contours, position = contour.legend.position, cex = cex.contours.legend)
  }

  fields::set.panel()
  
  # FINISCE FIGURA
  if(!is.null(file.name)) {dev.off(); message('figure saved in ', file.name)}

}


#' fa un dump del file .3d in una serie di immagini png corrispondenti a ciascuna
#' slice in z.
#'
#' @family display slices
#' @export
dumpslices.3d.2.png <- function(file.name,
                                x.min=-Inf, x.max=+Inf,
                                y.min=-Inf, y.max=+Inf,
                                z.min=-Inf, z.max=+Inf) {

  cat('reading values:', file.name, '\n')

  # leggi file 3d
  # apri connsessione
  file.3d <- file(file.name, "rb") # read binary

  # leggi l'header
  cat('reading header...\n')
  myline <- readLines(file.3d, n=8) # legge le prime 8 linee, ogni linea e' un elemento del vettore

  # parsing dell'header
  # splitta la stringa in sottostringhe delimitate da uno o piu' spazi (regexpr: " +")
  myline.splitted <- unlist(strsplit(myline[1], ' +'))
  Nx <- as.numeric(myline.splitted[1])
  myline.splitted <- unlist(strsplit(myline[2], ' +'))
  x <- as.numeric(myline.splitted)

  myline.splitted <- unlist(strsplit(myline[3], ' +'))
  Ny <- as.numeric(myline.splitted[1])
  myline.splitted <- unlist(strsplit(myline[4], ' +'))
  y <- as.numeric(myline.splitted)

  myline.splitted <- unlist(strsplit(myline[5], ' +'))
  Nz <- as.numeric(myline.splitted[1])
  myline.splitted <- unlist(strsplit(myline[6], ' +'))
  z <- as.numeric(myline.splitted)

  myline.splitted <- unlist(strsplit(myline[7], ' +'))
  Nv <- as.numeric(myline.splitted[1])
  values <- myline.splitted[2:length(myline.splitted)]

  Ntot <- Nx * Ny * Nz * Nv
  cat('number of voxels:', Nv, 'x', Nx, 'x', Ny, 'x', Nz, '=', Ntot, '\n')
  cat('variables:', values, '\n')

  # trasforma intervalli in coordinate puntuali
  x <- (x[1:Nx] + x[2:(Nx+1)])/2
  y <- (y[1:Ny] + y[2:(Ny+1)])/2
  z <- (z[1:Nz] + z[2:(Nz+1)])/2

  # crea e legge array (4d)
  cat('reading binary data...\n')
  Values.3d <- array(readBin(file.3d, numeric(), Ntot), dim=c(Nv, Nx, Ny, Nz))


  # selezione sottovolume

  cat('selecting subvolume...\n')

  # trova estremi nelle coordinate disponibili
  xx.min <- min(x[x>=x.min])
  yy.min <- min(y[y>=y.min])
  zz.min <- min(z[z>=z.min])
  xx.max <- max(x[x<=x.max])
  yy.max <- max(y[y<=y.max])
  zz.max <- max(z[z<=z.max])

  i.min <- which(x==xx.min)
  j.min <- which(y==yy.min)
  k.min <- which(z==zz.min)
  i.max <- which(x==xx.max)
  j.max <- which(y==yy.max)
  k.max <- which(z==zz.max)

  cat('subvolume: (', xx.min, ', ', xx.max, ') ',
      '(', yy.min, ', ', yy.max, ') ',
      '(', yy.min, ', ', yy.max, ')\n', sep='')
  cat('indexes: (', i.min, ', ', i.max, ') ',
      '(', j.min, ', ', j.max, ') ',
      '(', k.min, ', ', k.max, ')\n', sep='')

  # seleziona sottoarray
  Values.3d <- Values.3d[ , i.min:i.max, j.min:j.max, k.min:k.max]

  # trim delle coordinate
  x <- x[x>=xx.min & x<=xx.max]
  y <- y[y>=yy.min & y<=yy.max]
  z <- z[z>=zz.min & z<=zz.max]

  Nx <- length(x)
  Ny <- length(y)
  Nz <- length(z)


  # salva png nel folder

  cat('saving slices to pngs...\n')
  library(png)

  myfolder <- paste(file.name, '.slices.png', sep='')
  dir.create(myfolder)

  
  pb <- txtProgressBar(min = 0, max = Nv, style = 3)
  for (v in 1:Nv) {
    setTxtProgressBar(pb, v)
    #cat('saving pngs for variable:', values[v], '...\n')
    for (k in 1:Nz) {
      png.name <- paste(values[v], '_', z[k], 'y.png', sep='')
      # sanitize
      png.name <- gsub('/', '_', png.name)
      #cat (png.name, '\n')
      if (Nv > 1) {im.slice <- t(Values.3d[v, , , k])}
      else {
        #cat(dim(Values.3d), '\n')
        im.slice <- t(Values.3d[ , , k])
      }
      v.max <- max(im.slice)
      v.min <- min(im.slice)
      im.slice <- (im.slice - v.min)/(v.max - v.min)
      writePNG(im.slice, target=paste(myfolder, '/', png.name, sep=''))
    }
  }
  close(pb)
}


# REPORTING FUNCTIONS ----------------------------------------------------------


#' procedura per creare un display complessivo da un "plan"
display.all.plan <- function(plan)
{

  # recupera dati
  .ct <- get.ct(plan)
  .roi <- get.contours(plan)
  .values <- get.values(plan)
  .vois <- get.vois(plan)

  # display
  for(v in 1:.values$Nv) {
    display.slice.all(ct=.ct,
                      contours=.roi,
                      values=.values,
                      variable=.values$variables[v],
                      plan=plan)
  }
}


# DVH DISPLAY FUNCTIONS --------------------------------------------------------


# plot del dvh per una variabile specifica
# display.dvh.old <- function(dvh, plan=NULL, Diff=FALSE) {
# 
#   # check per vedere se è una lista
#   if (class(dvh[[1]]) != "numeric") {N <- length(dvh)} else {N <- 1}
# 
#   # text
#   if(N==1) {
#     if(!is.null(plan)) {main <- paste(plan$name, '-', dvh$voi)}
#     else {main <- dvh$voi}
#   } else {
#     if(!is.null(plan)) {main <- paste(plan$name, '- DVHs')}
#     else {main <- 'DVHs'}
#     #my.cols <- sample(colors(), N)
#     palette(rainbow(round(N*1.5))) # prende solo la prima parte del rainbow
#     my.cols <- palette()[1:N]
#   }
# 
#   # plot
#   if(N==1) {
#     if(!Diff) {
#       plot(dvh$value, dvh$volume, type='l', main=main, xlab=dvh$variable, ylab='%Volume')
#     } else {
#       hist(dvh$value, breaks=100, freq=FALSE, main=main, xlab=dvh$variable, ylab='Normalized Volume')
#     }
#   } else {
#     # check variabili
#     val.min <- 1e10
#     val.max <- -1e10
#     vois <- variables <- rep(' ', N)
#     for(i in 1:N) {
#       variables[i] <- dvh[[i]]$variable
#       vois[i] <- dvh[[i]]$voi
#       val.min <- min(val.min, min(dvh[[i]]$value))
#       val.max <- max(val.max, max(dvh[[i]]$value))
#     }
#     if(length(unique(variables))>1) {
#       message('error: inconsistent variables...')
#       return()
#     } else {message('dvhs variable: ', unique(variables))}
# 
#     # Add extra space to right of plot area; change clipping to figure
#     par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
# 
#     plot(dvh[[1]]$value, dvh[[1]]$volume,
#          type='l', col=my.cols[1],
#          xlim=c(val.min, val.max),
#          ylim=c(0, 1),
#          main=main, xlab=dvh[[1]]$variable, ylab='%Volume')
#     for(i in 2:N) {
#       lines(dvh[[i]]$value, dvh[[i]]$volume, col=my.cols[i])
#     }
#     legend("topright", inset=c(-0.5,0), legend=vois, col=my.cols, lty=1, title="VOIs", cex=0.6, bty='n')
#   }
# }


#' Display DVHs
#'
#' Display a single DVH or a list of DVHs
#'
#' @param dvh single DVH (created via \code{\link{dvh.evaluate}}) or a list of DVHs
#' @param plan plan object (optional)
#' @param Diff display "differential" DVH (boolean, optional)
#' @param alpha.color opacity of the plot (optional)
#' @param title the title of the plot (optional)
#' @param show.plot if \code{TRUE} it display the plot on screen. If \code{FALSE} it returns a ggplot plot structure if (optional)
#' @param decimate if the DVH points are too much (> 1000), it plots only a sample of them, to speed up the visualization (boolean, optional)
#' @param show.prescription plot the prescription over the DVH. The plan object needs to be specified also (boolean, optional)
#' @param filename the name of the file in which to save the figure (optional)
#' @param height,weight sizes (inches) of the figure to be saved in the file (optional)
#' @param fixed.scale Forces a common fixed scale when faceting (if there are different variables).
#' @param original.color use ggplot default colors
#' @param return.dataframe Return the data.frame of the DVH(s).
#'
#' @return If \code{show.plot} is \code{FALSE}, it returns a ggplot2 plot structure to be used for further processing. If \code{return.dataframe = TRUE} it returns the generated data.frame of the DVH(s).
#'
#' @family display dvh
#' @export
#' @import ggplot2
display.dvh <- function(dvh, plan=NULL,
                        Diff=FALSE,
                        alpha.color=1,
                        title=NULL,
                        show.plot=TRUE,
                        decimate=TRUE,
                        show.prescription=FALSE,
                        file.name=NULL,
                        height=7,
                        width=7,
                        original.colors=TRUE,
                        fixed.scale=FALSE,
                        return.dataframe=FALSE) {

  # usa le librerie ggplot2 (per fare prima...)
  #suppressMessages(library(ggplot2))

  # massimo numero di pti per dvh da visualizzare
  max.v <- 1000

  # check per vedere se è una lista
  if (class(dvh[[1]]) != "numeric") {N <- length(dvh)}
  else {N <- 1; dvh <- list(dvh)}

  # text
  if(!is.null(title)) {
    main <- title
  } else if(N==1) {
    if(!is.null(plan)) {main <- paste(plan$name, '-', dvh$voi)}
    else {main <- dvh$voi}
  } else {
    if(!is.null(plan)) {main <- paste(plan$name, '- DVHs')}
    else {main <- 'DVHs'}
  }

  # crea dataframe per ggplot
  init=1
  for(i in 1:N) {
    message('dataframing dvh #:', i)
    if(length(dvh[[i]]$value)==0) {warning('void dvh!'); init <- init+1; next}
    df.tmp <- data.frame(volume=dvh[[i]]$volume,
                         value=dvh[[i]]$value,
                         voi=dvh[[i]]$voi,
                         variable=dvh[[i]]$variable,
                         Nvoxel=dvh[[i]]$Nvoxel,
                         id=i)
    # riduce il numero di punti...
    if(decimate) {
      nr <- nrow(df.tmp)
      if(nr>max.v+2) {
        #message('decimating dvh...')
        index <- c(1, sort(sample(x=2:(nr-1), size=max.v)), nr)
        df.tmp <- df.tmp[index,]
      }
    }

    if(i==init) {df <- df.tmp} else (df <- rbind(df, df.tmp))
  }

  if(return.dataframe) {return(df)}

  #print(summary(df))
  unique.voi <- unique(df$voi)
  Ncol <- length(unique.voi)
  #my.cols <- sample(colors(), N)
  #palette(rainbow(round(Ncol*1.5))) # prende solo la prima parte del rainbow
  #my.cols <- palette()[1:Ncol]
  my.cols <- rainbow(round(Ncol*1.5))[1:Ncol]
  
  # check per vedere se ci sono display.color predefiniti.
  for(ic in 1:N) {
    if(!is.null(dvh[[ic]]$display.color)) {
      message('using display.color for ', dvh[[ic]]$voi)
      icol <- which(dvh[[ic]]$voi==unique.voi)
      if(dvh[[ic]]$display.color=="#FFFFFF"){dvh[[ic]]$display.color <- "#000000"}
      my.cols[icol] <- dvh[[ic]]$display.color
    }
  }
  

  # check per vedere quante variabili ci sono nella lista
  variables <- unique(df$variable)
  Nv <- length(variables)

  # scaling
  my.scale <- 1/max.v*100

  # prescrizione
  if(show.prescription) {
	  pres <- get.prescription(plan)[c(1,3,4,5,6)]
	  names(pres) <- c('voi',  'type', 'variable', 'value.pres', 'volumeFraction')
	  pres <- subset(pres, variable %in% unique(df$variable) & voi %in% unique(df$voi))
	  print(summary(pres))
  }

  # plot
  if(!Diff) {
    p <- ggplot(df) +
      geom_line(alpha=alpha.color, aes(x=value, y=volume*100, colour=voi, group=id)) +
      labs(y='%Volume', title=main, colour='VOI') +
      my.ggplot.theme()
      if(show.prescription) {
	      p <- p + geom_point(data=pres, aes(x=value.pres, y=volumeFraction*100, colour=voi, shape=type))
	}
  } else {
    alpha.color <- alpha.color/(N+1)
    p <- ggplot(df) +

      # solo Diff
      stat_density(position='dodge', alpha=alpha.color, aes(x=value, colour=voi, fill=voi, group=id)) +

      # both Diff+integral
      #stat_density(position='dodge', alpha=alpha.color, aes(x=value, y=..scaled..*25, colour=voi, fill=voi, group=id)) +
      #geom_line(aes(x=value, y=volume*100, colour=voi, group=id)) +

     labs(y='Normalized Volume', title=main, colour='VOI', fill='VOI') +
      scale_fill_manual(values=my.cols) +
      my.ggplot.theme() #+ coord_cartesian(xlim = c(1, 1.25))
  }
  
  if(!original.colors) {
    p <- p + scale_color_manual(values=my.cols) +
      scale_fill_manual(values=my.cols)
  }

  if(Nv==1) {
    p <- p + labs(x=variables)
  } else {
    if(fixed.scale) {
      p <- p + facet_wrap(~variable)
    } else {
      p <- p + facet_wrap(~variable, scales='free_x')
    }
  }

  if(!is.null(file.name)) {ggsave(plot=p, filename=file.name, height=height, width=width)}
  if(show.plot) {print(p)} else {return(p)}
}


#' Display DVHs
#'
#' Display a single DVH or a list of DVHs with the possibility to use cumulative or differential representations.
#'
#' @param dvh single DVH (created via \code{\link{dvh.evaluate}}) or a list of DVHs
#' @param type type of display ('cumulative', 'differential', or 'both')
#' @param alpha.color opacity of the plot (optional)
#' @param title the title of the plot (optional)
#' @param show.plot if \code{TRUE} it display the plot on screen. If \code{FALSE} it returns a ggplot plot structure if (optional)
#' @param decimate if the DVH points are too much (> 1000), it plots only a sample of them, to speed up the visualization (boolean, optional)
#' @param filename the name of the file in which to save the figure (optional)
#' @param height,weight sizes (inches) of the figure to be saved in the file (optional)
#' @return If \code{show.plot} is \code{FALSE}, it returns a ggplot2 plot structure to be used for further processing. If \code{return.dataframe = TRUE} it returns the generated data.frame of the DVH(s).
#'
#' @family display dvh
#' @export
#' @import ggplot2
display.dvh.combined <- function(dvh,
                                 type='cumulative', # scelta tra 'cumulative', 'differential', 'both'
                                 alpha.color=1,
                                 title=NULL,
                                 decimate=TRUE,
                                 unique.variable=FALSE,
                                 show.plot=TRUE,
                                 original.colors=TRUE,
                                 filename=NULL, width=7, height=7) {
  
  # usa le librerie ggplot2 (per fare prima...)
  #suppressMessages(library(ggplot2))
  
  # massimo numero di pti per dvh da visualizzare
  max.v <- 1000
  
  # check per vedere se è una lista
  if (class(dvh[[1]]) != "numeric") {N <- length(dvh)}
  else {N <- 1; dvh <- list(dvh)}
  
  # title
  # if(!is.null(title)) {
  #   main <- title
  # } else if(N==1) {
  #   if(!is.null(plan)) {main <- paste(plan$name, '-', dvh$voi)}
  #   else {main <- dvh$voi}
  # } else {
  #   if(!is.null(plan)) {main <- paste(plan$name, '- DVHs')}
  #   else {main <- 'DVHs'}
  # }
  main <- title
  
  # crea dataframe per ggplot
  init=1
  df.tmp <- df.tmp.diff <- NULL
  for(i in 1:N) {
    message('dataframing dvh #:', i)
    
    # cumulativo
    if(type=='cumulative' | type=='both') {
      if(length(dvh[[i]]$value)==0) {warning('void dvh!'); init <- init+1; next}
      df.tmp <- data.frame(volume=dvh[[i]]$volume*100,
                           value=dvh[[i]]$value,
                           voi=dvh[[i]]$voi,
                           variable=dvh[[i]]$variable,
                           Type='cumulative',
                           alpha=0,
                           id=i)
      # riduce il numero di punti...
      if(decimate) {
        nr <- nrow(df.tmp)
        if(nr>max.v+2) {
          #message('decimating dvh...')
          index <- c(1, sort(sample(x=2:(nr-1), size=max.v)), nr)
          df.tmp <- df.tmp[index,]
        }
      }
    }
    
    # differenziale
    if(type=='differential' | type=='both') {
      # stima della densità
      dens <- density(dvh[[i]]$value,
                      na.rm=TRUE,
                      bw='nrd0',
                      from=min(dvh[[i]]$value, na.rm=TRUE),
                      to=max(dvh[[i]]$value, na.rm=TRUE))
      df.tmp.diff <- data.frame(volume=dens$y,
                                value=dens$x,
                                voi=dvh[[i]]$voi,
                                variable=dvh[[i]]$variable,
                                Type='differential',
                                alpha=alpha.color/(N+1),
                                id=i)
    }
    
    if(i==init) {df <- rbind(df.tmp, df.tmp.diff)} else (df <- rbind(df, df.tmp, df.tmp.diff))
  }
  
  # check per vedere quante variabili ci sono nella lista
  variables <- unique(df$variable)
  Nv <- length(variables)
  
  # check per vedere se ci sono VOI ripetuti.
  if(unique.variable & length(unique(df$voi)) < N) {
    df$voi <- as.factor(paste(df$voi, df$variable))
  }
  #print(summary(df))
  #print(unique(df$id))
  
  # scala colori
  Ncol <- length(unique(df$voi))
  #my.cols <- sample(colors(), N)
  #palette(rainbow(round(Ncol*1.5))) # prende solo la prima parte del rainbow
  #my.cols <- palette()[1:Ncol]
  my.cols <- rainbow(round(Ncol*1.5))[1:Ncol]
  
  # plot
  p <- ggplot(df) +
    my.ggplot.theme() + #+ coord_cartesian(xlim = c(1, 1.25))
    geom_ribbon(aes(x=value, ymax=volume, ymin=0, fill=voi, alpha=Type, group=id)) +
    geom_line(aes(x=value, y=volume, colour=voi, group=id)) +
    labs(y='Normalized Volume', title=main, colour='VOI', fill='VOI')
  if(!original.colors) {
    p <- p + scale_color_manual(values=my.cols) +
      scale_fill_manual(values=my.cols)
  }
  p <- p + guides(alpha=FALSE, fill=FALSE)
  
  
  if(type=='both') {
    p <- p + scale_alpha_discrete(range=c(0, alpha.color/(N+1)))
    if(Nv==1) {
      p <- p + labs(x=variables) + facet_grid(Type~., scales='free_y')
    } else {
      if(!unique.variable) p <- p + facet_grid(Type~variable, scales='free')
      if(unique.variable) p <- p + labs(x=variables[1]) + facet_grid(Type~., scales='free_y')
    }
  }
  else {
    if(type=='differential') {p <- p + scale_alpha_discrete(range=c(alpha.color/(N+1), alpha.color/(N+1)))}
    else {p <- p + scale_alpha_discrete(range=c(0, 0))}
    if(Nv==1) {
      p <- p + labs(x=variables)
    } else {
      if(!unique.variable) p <- p + facet_grid(.~variable, scales='free')
      if(unique.variable) p <- p + labs(x=variables[1])
    }
  }
  
  if(!is.null(filename)) {ggsave(plot=p, filename=filename, width=width, height=height)}
  if(show.plot) {print(p)} else {return(p)}
}




#' plot di "dvh2D per un voi e due variabili specifiche
#'
#' @import ggplot2
#' @export
display.dvh2d <- function(values=values, vois=vois, variables=c('Dose[Gy]', 'Dose[Gy]'), voi=voi, alpha=0.2, means=FALSE, x.lim=NULL, y.lim=NULL) {

  #library(grid)
  #library(gtable)

  # crea d.f
  df1 <- dataframe.from.values(values=values, vois=vois, variables=variables[1], rois=voi)
  df2 <- dataframe.from.values(values=values, vois=vois, variables=variables[2], rois=voi)

  d.f <- data.frame(x=df1$value, y=df2$value)

  # Main scatterplot
  p1 <- ggplot(d.f, aes(x, y)) +
    geom_point(alpha=alpha) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    expand_limits(y = c(min(d.f$y) - 0.1 * diff(range(d.f$y)),
                        max(d.f$y) + 0.1 * diff(range(d.f$y)))) +
    expand_limits(x = c(min(d.f$x) - 0.1 * diff(range(d.f$x)),
                        max(d.f$x) + 0.1 * diff(range(d.f$x)))) +
    theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines")) +
    labs(x=variables[1], y=variables[2])
  if(means) {p1 <- p1 + geom_point(x=mean(d.f$x), y=mean(d.f$y), color='red')}
  if(!is.null(x.lim)) {p1 <- p1 + scale_x_continuous(limits=x.lim)}
  if(!is.null(y.lim)) {p1 <- p1 + scale_y_continuous(limits=y.lim)}

  # Horizontal marginal density plot - to appear at the top of the chart
  p2 <- ggplot(d.f, aes(x=x)) +
    geom_density(alpha=alpha, trim=TRUE, fill='black') +
    scale_x_continuous(expand = c(0, 0)) +
    expand_limits(x = c(min(d.f$x) - 0.1 * diff(range(d.f$x)),
                        max(d.f$x) + 0.1 * diff(range(d.f$x)))) +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = unit(c(1, 0.2, -0.5, 0.5), "lines"))
  if(means) {p2 <- p2 + geom_vline(xintercept=mean(d.f$x), color='red')}
  if(!is.null(x.lim)) {p2 <- p2 + scale_x_continuous(limits=x.lim)}

  # Vertical marginal density plot - to appear at the right of the chart
  p3 <- ggplot(d.f, aes(x=y)) +
    geom_density(alpha=alpha, trim=TRUE, fill='black') +
    scale_x_continuous(expand = c(0, 0)) +
    expand_limits(x = c(min(d.f$y) - 0.1 * diff(range(d.f$y)),
                        max(d.f$y) + 0.1 * diff(range(d.f$y)))) +
    coord_flip() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = unit(c(0.2, 1, 0.5, -0.5), "lines"))
  if(means) {p3 <- p3 + geom_vline(xintercept=mean(d.f$y), color='red')}
  if(!is.null(y.lim)) {p3 <- p3 + scale_x_continuous(limits=y.lim)}

  # Get the gtables
  gt1 <- ggplot_gtable(ggplot_build(p1))
  gt2 <- ggplot_gtable(ggplot_build(p2))
  gt3 <- ggplot_gtable(ggplot_build(p3))

  # Get maximum widths and heights for x-axis and y-axis title and text
  maxWidth <- unit.pmax(gt1$widths[2:3], gt2$widths[2:3])
  maxHeight <- unit.pmax(gt1$heights[4:5], gt3$heights[4:5])

  # Set the maximums in the gtables for gt1, gt2 and gt3
  gt1$widths[2:3] <- as.list(maxWidth)
  gt2$widths[2:3] <- as.list(maxWidth)

  gt1$heights[4:5] <- as.list(maxHeight)
  gt3$heights[4:5] <- as.list(maxHeight)

  # Combine the scatterplot with the two marginal boxplots
  # Create a new gtable
  gt <- gtable(widths = unit(c(7, 2), "null"), height = unit(c(2, 7), "null"))

  # Instert gt1, gt2 and gt3 into the new gtable
  gt <- gtable_add_grob(gt, gt1, 2, 1)
  gt <- gtable_add_grob(gt, gt2, 1, 1)
  gt <- gtable_add_grob(gt, gt3, 2, 2)

  # And render the plot
  grid.newpage()
  grid.draw(gt)
}

#' plot multiplo di dvh2D
#'
#' utilizza una lista di values, di vois, e di voi. C'è flessibilità su come
#' si può combinare la molteplicità: ad es. il nome del voi specificato
#' può essere comune a tutti i vois, oppure si può specificare per ciascun
#' values. E' necessario fornire un vettore per la legenda, per identificare
#' a mano i diversi contributi.
#'
#' E' possibile passare dei dataframe di parametri per ogni distribuzione
#' per i modelli da usare (a cui sono associati coppie di parametri specifiche):
#' - (Dose, Survival[.LM,.cMKM]) -> model.LQ
#' - (LETd, alpha) -> model.LM, model.cMKM, model.MKM
#' - i dataframe dei parametri devono avere la forma:
#' par1, par2, ... , id
#' dove il numero di righe = length(legend) e id[i] = legend[i]
#'
#' i parametri hanno i nomi:
#' model.LQ -> (alpha, beta)
#' model.LM -> (alpha0, m)
#' model.cMKM -> (alphaX, betaX, Rn, Rd)
#' model.MKM -> (alphaX, betaX, Rn, Rd)
#'
#' se è RBE.alpha=TRUE, viene anche usata l'alphaX (specificata nel modello) per calcolare l'RBE.alpha
#' del modello
#'
#' Nota: per il modello cMKM occorre specificare obbligatoriamente l'intervallo xlim
#'
#' @param legend.position The position of the legend. legend.position = 0 means no legend.
#'
#'
#' @import ggplot2 gtable grid
#' @export
display.dvh2d.multiple <- function(values=values, vois=vois, variables=c('Dose[Gy]', 'Dose[Gy]'), voi=voi,
				   alpha=0.2, means=FALSE, x.lim=NULL, y.lim=NULL,
				   model.LQ=NULL, model.LM=NULL, model.cMKM=NULL, model.MKM=NULL, RBE.alpha=FALSE,
				   legend=c('plan'), legend.position=1, different.model.colors=FALSE, file.name=NULL, height=7, width=7)
{

  #library(grid)
  #library(gtable)
  my.ggplot.theme(size=14)

  if(is.null(names(values))) {Nval <- length(values)} else {Nval <- 1}
  if(is.null(names(vois))) {Nvois <- length(vois)} else {Nvois <- 1}
  Nvoi <- length(voi)
  NN <- max(Nval, Nvois, Nvoi)
  alpha <- alpha/sqrt(NN)

  message('Nval:', Nval, ' Nvois:', Nvois, ' Nvoi:', Nvoi)

  # crea d.f
  for(i in 1:NN) {

    if(Nval>1) {my.values <- values[[i]]} else {my.values <- values}
    if(Nvois>1) {my.vois <- vois[[i]]} else {my.vois <- vois}
    if(Nvoi>1) {my.voi <- voi[i]} else {my.voi <- voi}

    df1 <- dataframe.from.values(values=my.values,
                                     vois=my.vois,
                                     variables=variables[1],
                                     rois=my.voi)
    df2 <- dataframe.from.values(values=my.values,
                                     vois=my.vois,
                                     variables=variables[2],
                                     rois=my.voi)
    df.tmp <- data.frame(x=df1$value, y=df2$value, id=legend[i])
    if(i==1) {d.f <- df.tmp} else {d.f <- rbind(d.f, df.tmp)}
  }

  #if(return.dataframe) {return(d.f)}

  # Main scatterplot
  p1 <- ggplot(d.f, aes(x, y, colour=as.factor(id))) +
    geom_point(alpha=alpha/2) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    expand_limits(y = c(min(d.f$y) - 0.1 * diff(range(d.f$y)),
                        max(d.f$y) + 0.1 * diff(range(d.f$y)))) +
    expand_limits(x = c(min(d.f$x) - 0.1 * diff(range(d.f$x)),
                        max(d.f$x) + 0.1 * diff(range(d.f$x)))) +
    theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines")) +
    labs(x=variables[1], y=variables[2], colour=NULL) +
    guides(colour=guide_legend(override.aes=list(alpha=1)))

    if(legend.position==1) {
      p1 <- p1 + theme(legend.justification=c(0,1), legend.position=c(0,1)) # in alto a sinistra
    } else if(legend.position==2) {
      p1 <- p1 + theme(legend.justification=c(1,1), legend.position=c(1,1)) # in alto a destra
    } else if(legend.position==3) {
      p1 <- p1 + theme(legend.justification=c(1,0), legend.position=c(1,0)) # in basso a destra
    } else if(legend.position==0) {
      p1 <- p1 + theme(legend.position="none")
    }

  # modello LM
  if(!is.null(model.LM)) {
    if(different.model.colors) {
      model.LM$id <- paste(model.LM$id, 'LM')
    }
    if(RBE.alpha) {
      model.LM$alpha0 <- model.LM$alpha0/model.LM$alpha.X
      model.LM$m <- model.LM$m/model.LM$alpha.X
    }
    p1 <- p1 + geom_abline(data=model.LM, aes(intercept=alpha0, slope=m, colour=id))
  }
  # modello cMKM
  if(!is.null(model.cMKM)) {
    letd.min <- x.lim[1]
    letd.max <- x.lim[2]
    if(different.model.colors) {
      model.cMKM$id <- paste(model.cMKM$id, 'cMKM')
    }
    for(iMKM in 1:nrow(model.cMKM)) {
      MKM.df.tmp <- alpha.beta.mkm(alphaX=model.cMKM$alphaX[iMKM],
                                   betaX=model.cMKM$betaX[iMKM],
                                   rN=model.cMKM$rN[iMKM],
                                   rd=model.cMKM$rd[iMKM],
                                   particleType=as.character(model.cMKM$particleType[iMKM]),
                                   lets=seq(letd.min, letd.max, length.out=30))
      MKM.df.tmp$id <- model.cMKM$id[iMKM]
      MKM.df.tmp$alpha.X <- model.cMKM$alpha.X[iMKM] # NOTA: OCCORRE RISISTEMRE UN PO' LE NOTAZIONI...
      if(iMKM==1) {MKM.df <- MKM.df.tmp} else {MKM.df <- rbind(MKM.df, MKM.df.tmp)}
    }
    if(RBE.alpha) {
      MKM.df$alpha <- MKM.df$alpha/MKM.df$alpha.X	# sostituisce per far prima
    }
    p1 <- p1 + geom_line(data=MKM.df, aes(x=let, y=alpha, colour=id, group=id))
  }
  if(means) {p1 <- p1 + geom_point(x=mean(d.f$x), y=mean(d.f$y), color='red')}
  if(!is.null(x.lim)) {p1 <- p1 + scale_x_continuous(limits=x.lim)}
  if(!is.null(y.lim)) {p1 <- p1 + scale_y_continuous(limits=y.lim)}
  #print(p1)

  # Horizontal marginal density plot - to appear at the top of the chart
  p2 <- ggplot(d.f, aes(x=x, colour=as.factor(id))) +
    geom_density(alpha=alpha, trim=TRUE, aes(fill=as.factor(id))) +
    scale_x_continuous(expand = c(0, 0)) +
    expand_limits(x = c(min(d.f$x) - 0.1 * diff(range(d.f$x)),
                        max(d.f$x) + 0.1 * diff(range(d.f$x)))) +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = unit(c(1, 0.2, -0.5, 0.5), "lines")) +
    scale_fill_discrete(guide=FALSE) +
    scale_colour_discrete(guide=FALSE)
  if(means) {p2 <- p2 + geom_vline(xintercept=mean(d.f$x), color='red')}
  if(!is.null(x.lim)) {p2 <- p2 + scale_x_continuous(limits=x.lim)}
  #print(p2)

  # Vertical marginal density plot - to appear at the right of the chart
  p3 <- ggplot(d.f, aes(x=y, colour=as.factor(id))) +
    geom_density(alpha=alpha, trim=TRUE, aes(fill=as.factor(id))) +
    scale_x_continuous(expand = c(0, 0)) +
    expand_limits(x = c(min(d.f$y) - 0.1 * diff(range(d.f$y)),
                        max(d.f$y) + 0.1 * diff(range(d.f$y)))) +
    coord_flip() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = unit(c(0.2, 1, 0.5, -0.5), "lines")) +
    scale_fill_discrete(guide=FALSE) +
    scale_colour_discrete(guide=FALSE)
  if(means) {p3 <- p3 + geom_vline(xintercept=mean(d.f$y), color='red')}
  if(!is.null(y.lim)) {p3 <- p3 + scale_x_continuous(limits=y.lim)}

  # Get the gtables
  gt1 <- ggplot_gtable(ggplot_build(p1))
  gt2 <- ggplot_gtable(ggplot_build(p2))
  gt3 <- ggplot_gtable(ggplot_build(p3))

  # Get maximum widths and heights for x-axis and y-axis title and text
  maxWidth <- unit.pmax(gt1$widths[2:3], gt2$widths[2:3])
  maxHeight <- unit.pmax(gt1$heights[4:5], gt3$heights[4:5])

  # Set the maximums in the gtables for gt1, gt2 and gt3
  gt1$widths[2:3] <- as.list(maxWidth)
  gt2$widths[2:3] <- as.list(maxWidth)

  gt1$heights[4:5] <- as.list(maxHeight)
  gt3$heights[4:5] <- as.list(maxHeight)

  # Combine the scatterplot with the two marginal boxplots
  # Create a new gtable
  gt <- gtable(widths = unit(c(7, 2), "null"), height = unit(c(2, 7), "null"))

  # Instert gt1, gt2 and gt3 into the new gtable
  gt <- gtable_add_grob(gt, gt1, 2, 1)
  gt <- gtable_add_grob(gt, gt2, 1, 1)
  gt <- gtable_add_grob(gt, gt3, 2, 2)

  # And render the plot
  if(!is.null(file.name)) {
    png(filename=file.name, width=width, height=height, res=300, units='in')
  }
  grid.newpage()
  grid.draw(gt)
  if(!is.null(file.name)) {
    dev.off()
  }
}


#' Visualizza una banda sul DVH
#'
#' accetta una lista "dvh.bands" (composta da due dvh: dvh.max e dvh.min)
#' se first.is.reference=TRUE allora visualizza il primo dvh della lista come riferimento (non è quindi usato per calcolare le bande)
#'
#' @export
#' @import ggplot2
display.dvh.bands <- function(dvh,
                              plan=NULL,
                              alpha=0.37,
                              eval.bands=TRUE,
                              alpha.color=0.25,
                              show.plot=TRUE,
                              file.name=NULL,
                              width=7,
                              height=7,
                              first.is.reference=FALSE,
                              with.mean=TRUE,
                              with.median=TRUE,
                              return.dataframe=FALSE) {

  # usa le librerie ggplot2 (per fare prima...)
  #suppressMessages(library(ggplot2))

  # se nella lista dvh ci sono più di 2 sv, allora calcola direttamente
  # le bande
  if(eval.bands) {
    message('evaluating bands...')
    if(first.is.reference) {
      dvh0 <- dvh[[1]]
      dvh <- dvh.evaluate.bands(dvh[2:length(dvh)], alpha)
    } else {
      dvh <- dvh.evaluate.bands(dvh, alpha)
    }
  }

  # crea dataframe per ggplot
  df <- data.frame(volume=dvh[[1]]$volume,
                   value.mean=dvh$dvh.mean$value,
                   value.median=dvh$dvh.median$value,
                   value.up=dvh$dvh.alpha.up$value,
                   value.lo=dvh$dvh.alpha.lo$value)

  if(!is.null(plan)) {
    tit <- paste(plan$name, ' - ', dvh[[1]]$voi, sep='')
  } else {tit <- dvh[[1]]$voi}

  p <- ggplot(df) +
    geom_ribbon(alpha=alpha.color, aes(x=volume*100, ymax=value.up, ymin=value.lo)) +
    coord_flip() +
    labs(title=tit, y=dvh[[1]]$variable, x='%Volume') +
    my.ggplot.theme()

  if(with.mean) {
    p <- p + geom_line(aes(x=volume*100, y=value.mean))
  }
  if(with.median) {
    p <- p + geom_line(linetype=2, aes(x=volume*100, y=value.median))
  }
  if(first.is.reference) {
    df0 <- data.frame(value=dvh0$value, volume=dvh0$volume)
    p <- p + geom_line(data=df0, colour='red', aes(x=volume*100, y=value))
  }

  if(!is.null(file.name)) {ggsave(plot=p, filename=file.name, width=width, height=height)}
  if(show.plot) {print(p)} else {return(p)}

}


#' Visualizza una banda sui DVH
#'
#' accetta una lista di "dvh"
#' i dvh vengono raggruppati a seconda del nome contenuto in dvh$voi.
#'
#' @export
#' @import ggplot2
display.dvh.bands.multiple <- function(dvh,
                                       dvh.reference=NULL,
                                       alpha=0.37,
                                       alpha.color=0.25,
                                       show.plot=TRUE,
                                       file.name=NULL,
                                       width=7,
                                       height=7,
                                       with.mean=TRUE,
                                       with.median=TRUE,
                                       title=NULL,
                                       return.dataframe=FALSE) {

  # usa le librerie ggplot2 (per fare prima...)
  #suppressMessages(library(ggplot2))

  # fa diventare il dvh.reference una lista
  if(!is.null(dvh.reference) & (!is.null(names(dvh.reference)) | length(names(dvh.reference)[1]=='value')!=0)) {
    message('one dvh for reference...')
    dvh.reference <- list(dvh.reference)
  }

  # identifica i nomi dei diversi voi:
  Ndvh <- length(dvh)
  voi <- rep('', Ndvh)
  for(i in 1:Ndvh) {
    voi[i] <- dvh[[i]]$voi
  }
  voi <- unique(voi)
  message('found vois: ', paste(voi, collapse=', '))


  # CREA dataframe DVH reference
  if(!is.null(dvh.reference)) {
    for(i in 1:length(dvh.reference)){
      df.reference.tmp <- data.frame(value=dvh.reference[[i]]$value,
                                     volume=dvh.reference[[i]]$volume,
                                     voi=dvh.reference[[i]]$voi)
      if(i==1) {df.reference <- df.reference.tmp} else {df.reference <- rbind(df.reference, df.reference.tmp)}
    }
  }

  # MAIN LOOP sui voi - DVHs
  for(vv in 1:length(voi)) {

    # identifica i DVH per il voi specifico
    dvh.tmp <-  list()
    index <- 1
    for(i in 1:Ndvh) {
      if(dvh[[i]]$voi==voi[vv]) {dvh.tmp[[index]] <- dvh[[i]]; index <- index+1}
    }

    # calcola direttamente
    # le bande
    message('evaluating bands...')
    dvh.bands <- dvh.evaluate.bands(dvh.tmp, alpha)

    # crea dataframe per ggplot
    df.tmp <- data.frame(volume=dvh.bands[[1]]$volume,
                         value.mean=dvh.bands$dvh.mean$value,
                         value.median=dvh.bands$dvh.median$value,
                         value.up=dvh.bands$dvh.alpha.up$value,
                         value.lo=dvh.bands$dvh.alpha.lo$value,
                         voi=voi[vv])

    if(vv==1) {df <- df.tmp} else {df <- rbind(df, df.tmp)}
  }
  
  if(return.dataframe) {return(list(df, df.reference))}

  p <- ggplot(df) +
    geom_ribbon(alpha=alpha.color, aes(x=volume*100, ymax=value.up, ymin=value.lo, fill=voi)) +
    coord_flip() +
    labs(title=title, y=dvh[[1]]$variable, x='%Volume') +
    my.ggplot.theme()


  if(with.mean) {
    p <- p + geom_line(aes(x=volume*100, y=value.mean, colour=voi))
  }
  if(with.median) {
    p <- p + geom_line(linetype=2, aes(x=volume*100, y=value.median, colour=voi))
  }
  if(!is.null(dvh.reference)) {
    p <- p + geom_line(data=df.reference, aes(x=volume*100, y=value, colour=voi))
  }


  if(!is.null(file.name)) {ggsave(plot=p, filename=file.name, width=width, height=height)}
  if(show.plot) {print(p)} else {return(p)}

}


# RENDERING 3D -----------------------------------------------------------------


#' rendering di isosuperfici da matrice values (openGL)
#'
#' @param mask a function of 3 arguments returning a logical array, a three dimensional logical array, or NULL. If not NULL, only cells for which mask is true at all eight vertices are used in forming the contour. Can also be a list of functions the same length as level.
#'
#' Rendering remoto usando Xvfb (da implementare)
#'
#' @export
# @import rgl misc3d
render.isosurfaces <- function(values, variable=NULL, levels=0, add=FALSE, alpha=NULL, color=NULL, file.name=NULL, axes=TRUE, mask=NULL, openGL=TRUE)
{
  
  # check per vedere se X11 è disponibile
  if(!capabilities(what='X11') | !openGL) {
    render.isosurfaces.static(values=values, variable=variable, levels=levels, add=add, alpha=alpha, color=color, file.name=file.name, axes=axes, mask=mask)
    return()
  }
  
  # carica esplicitamente le librerie
  library(rgl)
  library(misc3d)

  # identifica variabile
  v.index <- which(values$variable==variable)[1]
  Nv <- length(values$variables)
  if(Nv==1) {
    v.array <- values$values
  } else {
    v.array <- values$values[v.index,,,]
  }

  if(is.null(color)) {
    Ncol <- length(levels)
    my.cols <- rainbow(Ncol, start=0, end=0.6)[Ncol:1]
  } else {my.cols <- color}
  if(is.null(alpha)) {
    alpha <- 1/length(levels)
  }

  if(!is.null(mask)) {message('mask:');print(mask)}

  contour3d(v.array, levels, values$x, values$y, values$z, smooth=TRUE, add=add, alpha=alpha, color=my.cols, mask=mask)

  if(!add) {
    box3d()
    #rgl.viewpoint(30, 30, zoom=1, fov=30)
  }

  if(axes & !add) {axes3d(box=TRUE)}

  if(!is.null(file.name)) {snapshot3d(file.name)}

}

#' rendering di isosuperfici da matrice values
#'
#' @param mask a function of 3 arguments returning a logical array, a three dimensional logical array, or NULL. If not NULL, only cells for which mask is true at all eight vertices are used in forming the contour. Can also be a list of functions the same length as level.
#'
#' @export
# @import rgl misc3d
render.isosurfaces.static <- function(values, variable=NULL, levels=0, add=FALSE, alpha=NULL, color=NULL, file.name=NULL, axes=TRUE, mask=NULL)
{
  
  library(plot3D)
  
  #warning('render.isosurface.static not yet implemented')
  
  # identifica variabile
  v.index <- which(values$variable==variable)[1]
  Nv <- length(values$variables)
  if(Nv==1) {
    v.array <- values$values
  } else {
    v.array <- values$values[v.index,,,]
  }
  
  if(is.null(color)) {
    Ncol <- length(levels)
    my.cols <- rainbow(Ncol, start=0, end=0.6)[Ncol:1]
  } else {my.cols <- color}
  if(is.null(alpha)) {
    alpha <- 1/length(levels)
  }
  
  # contour3d(v.array, levels, values$x, values$y, values$z, smooth=TRUE, add=add, alpha=alpha, color=my.cols, mask=mask)
  isosurf3D(values$x, values$y, values$z, colvar=v.array, level = levels, col=my.cols, alpha=alpha, add=add)
  print(summary(v.array))
  print(add)
}


#' rendering 3D dei voi con isosuperfici
#'
#' @param mask a function of 3 arguments returning a logical array, a three dimensional logical array, or NULL. If not NULL, only cells for which mask is true at all eight vertices are used in forming the contour.
#'
#' @export

render.voi.isosurfaces <- function(vois=vois, voi=PTV, file.name=NULL, add=FALSE, alpha=NULL, mask=NULL, openGL=TRUE)
{

  vois.v <- vois

  Ncol <- length(voi)
  my.cols <- rainbow(round(Ncol*1.5))[1:Ncol]

  for(i in 1:length(voi))
  {
    message('evaluating isosurface for: ', voi[i])

    # crea la matrice values corrispondente al voi specificato per rendering con isosuperfici
    voi.index <- get.voi.logical(vois=vois, voi=voi[i])
    vois.v$values <- array(data=as.numeric(voi.index), dim=c(vois$Nx, vois$Ny, vois$Nz))

    vois.v$variables <- voi[i]

    if(i==1 & add==FALSE) {add <- FALSE} else {add <- TRUE}
    if(is.null(alpha)) {alpha <- 1/length(voi)}

    if(sum(vois.v$values)>0) {
      render.isosurfaces(values=vois.v, variable=voi[i], levels=0.5, add=add, alpha=alpha, color=my.cols[i], file.name=file.name, mask=mask, openGL=openGL)
    }
  }


}

#' mostra slice in maniera interattiva
#'
#' usa il pacchetto tkrplot (installa prima tk-dev con apt-get)
#' può visualizzare più di una variabile simultaneamente (mappandola sul "tempo")
#' (questo magari lo riformulerò quando userò immagini 4D...)
#'
#' @export
# @import tkrplot

display.slices.interactive <- function(values=values, variables=NULL, gray=FALSE)
{
 
  library(tkrplot) # carica esplicitamente la libreria
  
  # identifica variabili e formatta l'array temporale
  if(is.null(variables)) {variables <- values$variables}
  if(length(values$variables)==1) {
    v.array <- values$values
  } else if(length(variables)>1) {
    index.v <- which(values$variables==variables)
    v.array <- values$values[index.v,,,]
    v.array <- aperm(v.array, c(2,3,4,1))
  } else {
    index.v <- which(values$variables==variables)
    v.array <- values$values[index.v,,,]
  }

  if(!gray) {
    Ncol <- 256
    my.cols <- rainbow(Ncol, start=0, end=0.6)[Ncol:1]
  }

  slices3d(vol1=v.array, col1=my.cols)

}


# DISPLAY PROFILE FUNCTIONS ----------------------------------------------------

#' Display profile
#'
#' Deprecated: see "display.profiles"
#' @export
#' @family Profiles
display.profile <- function(profile.values,
                            profile.ct=NULL,
                            profile.names=NULL,
			                      depth.lim=NULL,
                            show.plot=TRUE,
                            file.name=NULL,
                            height=7,
                            width=7,
                            return.dataframe=FALSE) {

  #crea intervalli per la ct
  if(!is.null(profile.ct)) {
    names(profile.ct) <- c('variable', 'axis', 'depth', 'value.ct')
    d.depth <- mean(diff(profile.ct$depth))
    profile.ct$depth.min <- profile.ct$depth - d.depth/2
    profile.ct$depth.max <- profile.ct$depth + d.depth/2
    ct.variable <- unique(profile.ct$variable)
  }


  if(class(profile.values)=='list')
  {
    # combina eventuale lista profili in un unico dataframe
    for(i in 1:length(profile.values)) {
      profile.values.df.tmp <- profile.values[[i]]
      # strip away le colonne in eccesso (in maniera da creare profili consistenti)
      profile.values.df.tmp <- profile.values.df.tmp[c('variable', 'axis', 'depth', 'value')]
      if(!is.null(profile.names)) {
        profile.values.df.tmp$name <- profile.names[i]
      } else {
        profile.values.df.tmp$name <- as.factor(i)
      }
      profile.values.df.tmp$id <- i
      if(i==1) {profile.values.df <- profile.values.df.tmp} else {profile.values.df <- rbind(profile.values.df, profile.values.df.tmp)}
    }

    # estremi espliciti
    y.lim <- range(profile.values.df$value, na.rm = TRUE)
    if(diff(y.lim)==0) {y.lim <- c(y.lim[1]-0.5, y.lim[2]+0.5)}
    x.lim <- range(profile.values.df$depth, na.rm = TRUE)

    # variabile e asse del profilo
    values.variable <- unique(profile.values.df$variable)
    values.axis <- unique(profile.values.df$axis)

    #plot
    p <- ggplot()
    if(!is.null(profile.ct)) {
      p <- p + geom_rect(data=profile.ct, ymin=y.lim[1], ymax=y.lim[2], aes(xmin=depth.min, xmax=depth.max, fill=value.ct)) +
        scale_fill_gradient(low="gray", high='white') +
        labs(fill=ct.variable)
    }
    p <- p + geom_line(data=profile.values.df, aes(x=depth, y=value, colour=name, group=id)) +
      labs(x=values.axis, y=values.variable)
  }

  else

  {

    # estremi espliciti
    y.lim <- range(profile.values$value, na.rm = TRUE)
    if(diff(y.lim)==0) {y.lim <- c(y.lim[1]-0.5, y.lim[2]+0.5)}
    x.lim <- range(profile.values$depth, na.rm = TRUE)

    # variabile e asse del profilo
    values.variable <- unique(profile.values$variable)
    values.axis <- unique(profile.values$axis)


    p <- ggplot()
    if(!is.null(profile.ct)) {
      p <- p + geom_rect(data=profile.ct, ymin=y.lim[1], ymax=y.lim[2], aes(xmin=depth.min, xmax=depth.max, fill=value.ct)) +
        scale_fill_gradient(low="gray", high='white') +
        labs(fill=ct.variable)
    }
      p <- p + geom_line(data=profile.values, aes(x=depth, y=value)) +
      labs(x=values.axis, y=values.variable)
  }

  if(!is.null(depth.lim)) {p <- p + scale_x_continuous(limits=depth.lim)}

  my.ggplot.theme()
  if(!is.null(file.name)) {ggsave(plot=p, filename=file.name, height=height, width=width)}
  if(return.dataframe) {return(profile.values.df)}
  if(show.plot) {print(p)} else {return(p)}

}

#' Display profiles
#'
#' Display a single or multiple profiles. Optionally the profile could be plotted over a background color gradient representing the ct data (or other values).
#'
#' @param profile.values The main profiles dataframe to be plotted. It can contain a single profile or a multiple of profiles.
#' @param profile.ct the Background ct (or other data) profile.
#' @param show.ct.legend Boolean. A colorbar for the levels containted in profile.ct is displayed if TRUE.
#' @param x.axis.name Optional label for the x axis (default uses the data stored in profile.values if unique).
#' @param y.axis.name Optional label for the y axis (default uses the data stored in profile.values if unique).
#' @param name Optional label for the set of profiles.
#' @param depth.lim Optional limits for the x axis.
#' @param show.plot Show plot.
#' @param file.name File name for saving the plot.
#' @param height The height of the saved plot image (inches).
#' @param width The width of the saved plot image (inches).
#' @param return.plot return a ggplot2 plot if TRUE.
#' @export
#' @family Profiles
display.profiles <- function(profile.values,
                             profile.ct=NULL, show.ct.legend=TRUE,
                             y.axis.name=NULL,
                             x.axis.name=NULL,
                             name=NULL,
                             depth.lim=NULL,
                             file.name=NULL, height=7, width=7, dpi=300,
                             show.plot=TRUE, return.plot=FALSE)
{
  
  #crea intervalli per la ct
  if (!is.null(profile.ct)) {
    names(profile.ct) <- c("variable", "axis", "depth", "value.ct")
    d.depth <- mean(diff(profile.ct$depth))
    profile.ct$depth.min <- profile.ct$depth - d.depth/2
    profile.ct$depth.max <- profile.ct$depth + d.depth/2
    ct.variable <- unique(profile.ct$variable)[1]
    
    # estremi espliciti
    y.lim <- range(profile.values$value, na.rm = TRUE)
    if(diff(y.lim)==0) {y.lim <- c(y.lim[1]-0.5, y.lim[2]+0.5)} # nel caso di profilo "piatto"
    #x.lim <- range(profile.values$depth, na.rm = TRUE)
  }
  
  # variabile e asse del profilo
  if(is.null(x.axis.name)) {
    x.axis.name <- unique(profile.values$axis)
    if(length(x.axis.name)>1) x.axis.name <- 'depth'
  }
  variables <- unique(profile.values$variable)
  if(is.null(y.axis.name)) {
    y.axis.name <- variables
    if(length(y.axis.name)>1) y.axis.name <- 'value'
  }
  
  p <- ggplot()
  if(!is.null(profile.ct)) {
    p <- p + geom_rect(data=profile.ct, ymin=y.lim[1], ymax=y.lim[2], aes(xmin=depth.min, xmax=depth.max, fill=value.ct))
    if(show.ct.legend) p <- p + scale_fill_gradient(low="gray", high='white') + labs(fill=ct.variable)
    else p <- p + scale_fill_gradient(low="gray", high='white', guide=FALSE)
  }
  p <- p + geom_line(data=profile.values, aes(x=depth, y=value, colour=variable)) +
    labs(x=x.axis.name, y=y.axis.name, colour=name)
  if(length(variables)==1) p <- p + scale_colour_discrete(guide=FALSE)
  if(!is.null(depth.lim)) p <- p + scale_x_continuous(limits=depth.lim)
  
  my.ggplot.theme()
  if(!is.null(file.name)) {ggsave(plot=p, filename=file.name, height=height, width=width, dpi=dpi)}
  if(show.plot) {print(p)}
  if(return.plot) {return(p)}
}



# UTILITIES ====================================================================

#' Generate gray colormap for Hounsfield numbers
#'
#' @param HU.range The range c(HU.min, HU.max) of the Hounsfield numbers of the CT
#' @param HU.windos The range c(HU.min, HU.max) to map to c(0,1) gray map.
#' @export
colormap.ct <- function(HU.range=c(-1000, 3000), HU.window=c(-1000,3000))
{

  # controlla se la finestra è fuori range...
  scaling.HU.window <- FALSE
  if(HU.window[1]<HU.range[1]) {HU.window[1] <- HU.range[1]; scaling.HU.window <- TRUE}
  if(HU.window[2]>HU.range[2]) {HU.window[2] <- HU.range[2]; scaling.HU.window <- TRUE}
  if(scaling.HU.window) {message('HU.window scaled to HU.range')}



  Nc <- max(HU.range[2] - HU.range[1], 2)

  interval <- ((1:Nc)-1)/Nc

  if(diff(HU.range)==0) {
    col.ct <- 0.5
    return(col.ct)
  }


  wmin <- max(HU.window[1] - HU.range[1], 1)
  wmax <- min(HU.window[2] - HU.range[1], Nc)

  #message(wmin, ' ', wmax)
  #print(interval)

  interval[1:wmin] <- 0
  interval[(wmin+1):(wmax)] <- seq(0, 1, length.out=(wmax-wmin))
  interval[(wmax+1):Nc] <- 1

  col.ct <- gray(interval)

  return(col.ct)
}


#' imposta stile generale per ggplot2
#' @export
my.ggplot.theme <- function(size=14)
{
  # assume che le librerie ggplot siano già caricate
  my.theme <- theme_bw(size) + theme(axis.text = element_text(colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()); theme_set(my.theme)
  theme_set(my.theme)
  
  #pt <- theme(legend.key = element_rect(colour = 'white')) +
  #  theme(panel.border = element_rect(colour = "black"))

    # altre possibilità..
    # theme(panel.grid.major = element_line(colour = rgb(0.8, 0.8, 0.8)))
  #return(pt)
}


# BEAMS DISPLAY FUNCTIONS ------------------------------------------------------

#' Beam statistics plot
#'
#' @param beams The beams dataframe.
#' @param plan The plan object.
#' @param numeric Display a numeric ID for the fields (default uses a beamLine+angles code).
#' @param shot.plot Show plot.
#' @param file.name File name for saving the plot.
#' @param height The height of the saved plot image (inches).
#' @param width The width of the saved plot image (inches).
#' @param dpi dpi of the saved image.
#' @param field.names Optional vector of the names for the different fields
#' @return If show.plot=FALSE, it returns a ggplot object.
#' @export
#' @import ggplot2
#' @family Beams
display.beams <- function(beams,
                          plan=NULL,
                          numeric=FALSE,
                          show.plot=TRUE,
                          file.name=NULL,
                          height=7,
                          width=7,
                          dpi=300,
                          field.names=NULL)
{
  if(!is.null(plan)) {
    my.title <- paste(plan$name, 'Beams', sep=' - ')
  } else {
    my.title <- 'Beams'
  }

  # aggiunge field ID
  if(!('field' %in% colnames(beams))) {beams <- add.field(beams, numeric=numeric)}
  if(!is.null(field.names)) {beams$field <- field.names[beams$field]}

  # digits
  beams$energy.f <- as.factor(round(beams$energy, digits=1))

  # plot
  my.ggplot.theme()
  p <- ggplot(beams) +
    stat_count(aes(x=energy.f, weight=fluence, fill=field)) +
    coord_flip() +
    labs(y='Total number of primary ions', x='Energy Layers [MeV/u]', title=my.title) +
    theme(axis.title.y=element_text(vjust=.2))


  if(!is.null(file.name)) {ggsave(plot=p, filename=file.name, height=height, width=width, dpi = dpi)}
  if(show.plot) {print(p)} else {return(p)}
}


#' Display spots (3D)
#'
#' @param beams the beams dataframe
#' @param vois the vois object
#' @param voi(s) names
#' @param display.iso display the isocenter
#' @param alpha.spot the opacity value for the spot points
#'
#' @export
# @import rgl misc3d
display.spots <- function(beams, vois=NULL, voi=NULL, display.iso=TRUE, alpha.spot=0.2, alpha.voi=1)
{
  # plot 3d
  plot3d(beams$x_s, beams$y_s, beams$z_s, xlab='x [mm]', ylab='y [mm]', zlab='z [mm]', type='p', size=1, alpha.spot=alpha.spot)

  # isocentro
  if(display.iso) {
    iso <- aggregate(data.frame(N=rep(1,nrow(beams))), by=list(x=beams$x_iso, y=beams$y_iso, z=beams$z_iso), sum)
    print(iso)
    plot3d(iso$x, iso$y, iso$z, type='p', size=10, col='red', add=TRUE)
  }

  # vois
  if(!is.null(vois) & !is.null(voi)) {
    render.voi.isosurfaces(vois=vois, voi=voi, add=TRUE, alpha=alpha.voi)
  }
}


#' Display rays (3D)
#'
#' @param rays the rays dataframe.
#' @param alpha the opacity value for the point and lines.
#' @param ray.length Length of the rays. If ray.lengt=0 ti will plot only the spots.
#' @param add Add to existing 3D plot.
#' @export
# @import rgl misc3d
display.rays <- function(rays, alpha=1, ray.length=1, add=FALSE)
{
  # plot 3d
  plot3d(rays$X, rays$Y, rays$Z, xlab='x [mm]', ylab='y [mm]', zlab='z [mm]', type='p', size=2, alpha.spot=alpha, add=add, col='red')
  if(ray.length!=0)
  for(i in 1:nrow(rays)) {
    xx <- c(rays$X[i], rays$X[i] + rays$xn[i]*ray.length)
    yy <- c(rays$Y[i], rays$Y[i] + rays$yn[i]*ray.length)
    zz <- c(rays$Z[i], rays$Z[i] + rays$zn[i]*ray.length)
    plot3d(xx, yy, zz, type='l', add=TRUE, alpha=alpha/2)
  }
}


#' Beam-port splot
#'
#' @param beams beams dataframe.
#' @param plan The plan object.
#' @param numeric Display a numeric ID for the fields (default uses a beamLine+angles code).
#' @param shot.plot Show plot.
#' @param file.name File name for saving the plot.
#' @param height The height of the saved plot image (inches).
#' @param width The width of the saved plot image (inches).
#' @param dpi The dpi of the saved image.
#' @param field.names Optional vector of the names for the different fields
#' @return If show.plot=FALSE, it returns a ggplot object.
#' @export
#' @import ggplot2
#' @family Beams
display.beamports <- function(beams,
                          plan=NULL,
                          numeric=FALSE,
                          show.plot=TRUE,
                          file.name=NULL,
                          height=7,
                          width=7,
                          dpi=300,
                          field.names=NULL)
{
  if(!is.null(plan)) {
    my.title <- paste(plan$name, 'Beam-ports', sep=' - ')
  } else {
    my.title <- 'Beam-ports (fields)'
  }

  # aggiunge field ID
  beams <- add.field(beams, numeric=numeric)
  if(!is.null(field.names)) {beams$field <- field.names[beams$field]}

  beams.a <- aggregate(list(Npart=beams$fluence), list(deflX=beams$deflX, deflY=beams$deflY, beam.port=beams$field), sum)

  # plot
  p <- ggplot(beams.a) +
    geom_point(aes(x=deflX, y=deflY, colour=Npart, size=Npart)) +
    labs(y='deflY [mm]', x='deflX [mm]', coulour='N part.', title=my.title) +
    facet_wrap(~beam.port)
    #theme(axis.title.y=element_text(vjust=.2))
  my.ggplot.theme()

  if(!is.null(file.name)) {ggsave(plot=p, filename=file.name, height=height, width=width)}
  if(show.plot) {print(p)} else {return(p)}
}


# CONTOURS DISPLAY FUNCTIONS ---------------------------------------------------

#' Add a legend with contours names/colours.
#'
#' If colours are not specitied in the contours dataframe, they will be added using the "raimbow" colour selection.
#' @param contours The contours dataframe.
#' @param add Put the legend in the active plot. If FALSE, it returns a stand alone legend in a new display.
#' @param position The position of the legend ('topleft', 'bottomleft', etc.)
#' @param cex dimension of the text font
#' @param ... pther parameters to legend()
#' @export
#' @family Contours
display.contours.legend <- function(contours,
                                    add=TRUE,
                                    position='topleft',
                                    cex=0.8,
                                    ...)
{
  if( !('display.color' %in% colnames(contours)) ) {
    contours <- add.colours.contours(contours)
  }
  if(!add) {
    plot.new()
    position <- 'center'
  }
  legend(position, legend = unique(contours$contour), col = unique(contours$display.color), lty=1, cex=cex, ...)
}