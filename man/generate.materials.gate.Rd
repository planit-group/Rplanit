% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gate-functions.R
\name{generate.materials.gate}
\alias{generate.materials.gate}
\title{Generate Materials (Gate)}
\usage{

  generate.materials.gate(MaterialTable = "./data/Schneider2000MaterialsTable.txt",
  DensityTable = "./data/Schneider2000DensitiesTable.txt",
  DensityTolerance = 0.1,
  OutputMaterialDatabase = "./data/ct-HUmaterials.db",
  OutputHUMaterial = "./data/ct-HU2mat.txt")
}
\arguments{
\item{MaterialTable}{name of the file containing the conversion of CT number to elemental weight.}

\item{DensityTable}{name of the file containing the conversion of CT number to density values.}

\item{DensityTolerance}{the density step for the generation of the material database.}

\item{OutputMaterialDatabase}{output file name for the material database.}

\item{OutputHUMaterial}{output file name for the CT number to material conversion.}
}
\description{
Generate the material database for Gate simulations.
}
\seealso{
Other GateMaterials: \code{\link{read.materials.gate}}
}
\concept{GateMaterials}
