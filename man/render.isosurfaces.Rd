% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/display-functions.R
\name{render.isosurfaces}
\alias{render.isosurfaces}
\title{rendering di isosuperfici da matrice values (openGL)}
\usage{
render.isosurfaces(values, variable = NULL, levels = 0, add = FALSE,
  alpha = NULL, color = NULL, file.name = NULL, axes = TRUE,
  mask = NULL, openGL = TRUE)
}
\arguments{
\item{mask}{a function of 3 arguments returning a logical array, a three dimensional logical array, or NULL. If not NULL, only cells for which mask is true at all eight vertices are used in forming the contour. Can also be a list of functions the same length as level.

Rendering remoto usando Xvfb (da implementare)}
}
\description{
rendering di isosuperfici da matrice values (openGL)
}
