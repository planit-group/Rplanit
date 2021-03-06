% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/values-functions.R
\name{get.profiles}
\alias{get.profiles}
\title{Get profiles}
\usage{
get.profiles(values, variables = NULL, ...)
}
\arguments{
\item{values}{the \code{values object}}

\item{variables}{the variables to be profiled. If it is not specified, al the available variables in values will be used.}

\item{x, y, z}{definition of the axis (only two coordinates should be defined). For example \code{x=4.4} and \code{z=0.5} specifies a profile along the y-axis passing through x=4.4 and z=0.5}

\item{integrate}{perform an integration over the two coordinates (boolean, optional). Note: integration is not yet implemented around an arbitrary ray.}

\item{ray}{The ray object, Note the ray is a data.frame defined by 6 components (X, Y, Z, xn, yn, zn), i.e. point coordinates + normalized vector direction.}

\item{return.voxel.index}{Returns also the index of the voxel crossed by the ray (used only for a ray profile).}

\item{return.xyz}{Returns also the x,y,z of the points along the ray (used only for a ray profile).}
}
\value{
a dataframe containing the specified profile.
}
\description{
Get multiple profiles (values along a specific axis, or voxel line) from a \code{values} object. The axis is specifed trough the definition of two of the following three coordinates: x, y, z, or with a single ray object for arbitrary directions.
}
\concept{Values, Values Manipulation, Profiles}
