\name{display.profile}
\alias{display.profile}
\title{Display profile}
\usage{
display.profile(profile.values, profile.ct = NULL, profile.names = NULL,
  show.plot = TRUE, file.name = NULL, height = 7, width = 7)
}
\arguments{
  \item{profile.values}{the main profile dataframe to be
  plotted. It can be a single profile or a list of
  profiles}

  \item{profile.ct}{the background ct (or other data)
  profile}

  \item{profile.name}{a vector containing the names of the
  profiles. It is used if a list of different profiles is
  used as imput.}

  \item{shot.plot}{Show plot.}

  \item{file.name}{File name for saving the plot.}

  \item{height}{The height of the saved plot image
  (inches).}

  \item{width}{The width of the saved plot image (inches).}
}
\value{
If sho.plot=FALSE, it returns a ggplot object.
}
\description{
Display the profile. Optionally the profile could be
plotted over a background color gradient representing the
ct data (or other values).
}
\seealso{
Other Profiles: \code{\link{identify.slice}}
}

