\name{display.profile}
\alias{display.profile}
\title{Display profile}
\usage{
display.profile(profile.values, profile.ct = NULL, profile.names = NULL,
  depth.lim = NULL, show.plot = TRUE, file.name = NULL, height = 7,
  width = 7, return.dataframe = FALSE)
}
\arguments{
  \item{profile.values}{the main profile dataframe to be
  plotted. It can be a single profile or a list of
  profiles}

  \item{profile.ct}{the background ct (or other data)
  profile}

  \item{profile.names}{a vector containing the names of the
  profiles. It is used if a list of different profiles is
  used as input.}

  \item{depth.lim}{The x-axis range to display. (automatic
  if not specified).}

  \item{shot.plot}{Show plot.}

  \item{file.name}{File name for saving the plot.}

  \item{height}{The height of the saved plot image
  (inches).}

  \item{width}{The width of the saved plot image (inches).}

  \item{return.dataframe}{return the dataframe created for
  the plot (with ggplot2).}
}
\value{
If show.plot=FALSE, it returns a ggplot object. If
return.dataframe=TRUE, it returns the the dataframe created
for the plot (with ggplot2).
}
\description{
Display the profile. Optionally the profile could be
plotted over a background color gradient representing the
ct data (or other values).
}

