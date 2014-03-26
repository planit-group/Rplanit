# REPORTING FUNCTIONS ----------------------------------------------------------

#' Generate report
#' 
#' Generate a report summarizing the plan in a file with the markdown format (it can be converted in html via rstudio or the "markdown" package).
#' 
#' @param plan the \code{plan} object
#' @param file.name the file name
#' @param N.slice number of slices to plot for CT and value distributions
#' @param html output html file. Note: there is an issue with markdown package that prevents the correct latex equation display in html. Opening the markdown (*.md) file in Rstudio, and then previewing it, overcomes this issue (Rstudio uses an internal markdown to html conversion).
#' 
#' @export
#' @import knitr
generate.report <- function(plan, file.name='report', N.slice=6, html=FALSE)
{
  
  plan.name <- plan$name
  my.file.name <- 'report_temp'
  
text.rmd <- "
```{r echo=FALSE, results='hide', message=FALSE}
library(dektools)
plan <- read.plan('@plan.name@')
  
ct <- get.ct(plan)
contours <- get.contours(plan)
vois <- get.vois(plan)
values <- get.values(plan)
beams <- get.beams(plan)
  
N.slice <- @N.slice@
N.variables <- length(values$variables)
z.min <- min(contours$z, na.rm=TRUE); z.max <- max(contours$z, na.rm=TRUE)

ct.Dx <- diff(range(ct$x)); ct.Dy <- diff(range(ct$y))
d <- sqrt(36 / (ct.Dx*ct.Dy))
ct.fig.x <- d * ct.Dx * 1.15; ct.fig.y <- d * ct.Dy

values.Dx <- diff(range(values$x));values.Dy <- diff(range(values$y))
d <- sqrt(36 / (values.Dx*values.Dy))
values.fig.x <- d * values.Dx * 1; values.fig.y <- d * values.Dy
```
  
Plan: `r plan$name`
========================================================
  
  
CT
--------------------------------------------------------
CT information:
- CT size: `r max(ct$x)-min(ct$x)` $\\times$ `r max(ct$y)-min(ct$y)` $\\times$ `r max(ct$z)-min(ct$z)` mm$latex ^3$
- number of voxels: `r ct$Nx` $\\times$ `r ct$Ny` $\\times$ `r ct$Nz` $=$ `r ct$Nx*ct$Ny*ct$Nz`
- voxels size: `r mean(diff(ct$x))` $\\times$ `r mean(diff(ct$y))` $\\times$ `r mean(diff(ct$z))` mm$^3$
  
```{r fig.width=8, fig.height=6, echo=FALSE}
hist(ct$values, breaks=100, freq=FALSE, xlab=ct$variable, main='Histogram of Hounsfield Number')
```
  
```{r fig.width=ct.fig.x, fig.height=ct.fig.y, echo=FALSE, results='hide', message=FALSE}
z.slice <- seq(z.min, z.max, length.out=N.slice+2)[2:(N.slice+1)]
for(i in 1:N.slice) {
display.slice.ct(ct, roi=contours, z=z.slice[i])
}
```
  
Contours
--------------------------------------------------------
Number of VOIs: `r length(unique(contours$contour))`
  
```{r echo=FALSE, results='asis'}
contours.vol <- volume.contours(contours)
names(contours.vol)[2] <- 'volume[mm$^3$]'
contours.uni <- unique(contours[c('id', 'contour', 'tissue', 'type')])
contours.uni.vol <- merge(contours.uni, contours.vol)
kable(contours.uni.vol, format='markdown')
```

Fields
--------------------------------------------------------
Fields (beam-ports) information:
- Beamline model: `r plan$beamLine`
- Number of fields (beam-ports): `r nrow(plan$fields)`

```{r echo=FALSE, results='asis'}
fields.table <- plan$fields[c('targetVOI', 'iecGantryAngle', 'iecPatientSupportAngle', 'interSpotSpacing.x', 'interSpotSpacing.y', 'interSpotSpacing.z', 'spotsExtensionOutsideTarget')]
names(fields.table) <- c('target VOI', 'Gantry Angle [deg]', 'Patient Support Angle [deg]', 'Spot Spacing x', 'Spot Spacing y', 'Spot Spacing z', 'Spots Extension [mm]')
kable(fields.table, format='markdown')
```
(Note: spot spacing value is relative to the beam spot standard deviation.)

Fields (beam-ports) statistics:
- Number of particles: `r sum(beams$fluence)`
- Used energies: `r length(unique(beams$energy))` (from `r min(beams$energy)` to `r max(beams$energy)` MeV/u)

```{r fig.width=8, fig.height=(2+0.15*length(unique(beams$energy))), echo=FALSE, results='hide', message=FALSE}
display.beams(beams)
```

  
Values
--------------------------------------------------------
Values information:
- number of evaluated variables: `r length(unique(values$variables))`
- variable names: ` `r unique(values$variables)` `
- number of voxels of the computig grid: `r values$Nx` $\\times$ `r values$Ny` $\\times$ `r values$Nz` $=$ `r values$Nx*values$Ny*values$Nz`
- voxels size: `r mean(diff(values$x))` $\\times$ `r mean(diff(values$y))` $\\times$ `r mean(diff(values$z))` mm$^3$
  
```{r fig.width=values.fig.x, fig.height=values.fig.y, echo=FALSE, results='hide', message=FALSE, warning=FALSE}
z.slice <- seq(z.min, z.max, length.out=(N.slice+2))[2:(N.slice+1)]
for(j in 1:N.variables) {
  for(i in 1:N.slice) {
    display.slice.all(ct=ct,
                      contours=contours,
                      values=values,
                      variable=values$variables[j],
                      cont=FALSE,
                      z=z.slice[i])
  } 
}
```
  
DVHs
--------------------------------------------------------
```{r fig.width=(4+4*N.variables), fig.height=8, echo=FALSE, results='hide', message=FALSE}
for(i in 1:N.variables) {
  dvhs.tmp <- dvh.evaluate.all(values, vois, variable=values$variables[i])
  if(i==1) {dvhs <- dvhs.tmp} else {dvhs <- c(dvhs, dvhs.tmp)}
}
display.dvh.combined(dvhs)
```

Prescription
---------------------------------------------------------
Number of constraints: `r nrow(plan$prescription)`

```{r echo=FALSE, results='asis', message=FALSE, warning=FALSE}
prescription.table <- check.prescription(values=values, vois=vois, prescription=plan$prescription)
prescription.table <- prescription.table[, !'VOIIndex'==names(prescription.table)]
kable(prescription.table, format='markdown')
```
  
  
"
  
  # sostituzione con i parametri esterni
  text.rmd <- gsub(pattern='@plan.name@', replacement=plan.name, x=text.rmd, fixed=TRUE)
  text.rmd <- gsub(pattern='@N.slice@', replacement=N.slice, x=text.rmd, fixed=TRUE)
  #text.rmd <- gsub(pattern="\\", replacement='\\\\', x=text.rmd, fixed=TRUE)
  
  # file names
  my.file.name.rmd  <- paste(my.file.name, 'Rmd', sep='.')
  my.file.name.md  <- paste(my.file.name, 'md', sep='.')
  file.name.md <- paste(file.name, 'md', sep='.')
  #file.name.html <- paste(file.name, 'html', sep='.')
  
  # compile
  message('writing ', my.file.name.rmd)
  cat(text.rmd, file=my.file.name.rmd)
  #text.html <- knit2html(text=text.rmd)
  message('knitting ', my.file.name.rmd)
  knit(my.file.name.rmd, encoding='UTF-8')
  
  
  # converti in html
  if(html) {
    markdownToHTML(file=my.file.name.md)
  }
  
  
  # (re)move files
  file.rename(from=my.file.name.md, to=file.name.md)
  file.remove(my.file.name.rmd)
  
  #cat(text.html, file=file.name.html)
  
  #if(interactive()) browseURL(file.name.html)
  
  #return(text.rmd)
  
}