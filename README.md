# Rplanit

Radiotherapy with ion beams: planning, simulations and analysis

*Author*: Andrea Attili

*Maintainer*: Andrea Attili <andrea.attili@roma3.infn.it>

*Description*: A collection of methods and procedures to perform simulations and
analysis related to radiotherapy treatment planning (with ion beams). It includes a set
of radiobiological evaluations and display possibilities.


## Installation

NOTE: the following instructions are mostly outdated.

To perform radiobiological evaluations see also [Survival](https://github.com/batuff/Survival). An integration of Rplanit with Survival is underway.

To perform simulations and optimizations of treatments (TPS) an integration with the updated code Planit (based on [Russo et al. 2016](http://dx.doi.org/10.1088/0031-9155/61/1/183) and with [Topas](http://www.topasmc.org/) is underway.

### Ubuntu (12.04 64bit)
* The following libraries need to be installed from the Ubuntu repository (from the terminal shell prompt):
```
$ sudo apt-get install libcurl4-openssl-dev
$ sudo apt-get install tk8.5-dev
$ sudo apt-get install libsdl1.2-dev
$ sudo apt-get install mesa-common-dev # optional, for OpenGL rendering
$ sudo apt-get install libssl-dev
$ sudo ln -s /lib/x86_64-linux-gnu/libcrypto.so.1.0.0 /lib/x86_64-linux-gnu/libcrypto.so.6
$ sudo ln -s /lib/x86_64-linux-gnu/libssl.so.1.0.0 /lib/x86_64-linux-gnu/libssl.so.6
$ sudo apt-get install libgsl0-dev
$ sudo apt-get install liblog4cxx10-dev
```
* Install the latest version of R (for Ubuntu 12.04) following the instruction at this [link](http://livesoncoffee.wordpress.com/2012/12/09/installing-r-on-ubuntu-12-04/).
* Optional (but recommended): install the latest version of RStudio from this [link](http://www.rstudio.com/products/rstudio/download/).

### Ubuntu (14.04+ 64bit)
* The following libraries need to be installed from the Ubuntu repository (from the terminal shell prompt):
```
$ sudo apt-get install libcurl4-openssl-dev
$ sudo apt-get install tk8.6-dev
$ sudo apt-get install libsdl1.2-dev
$ sudo apt-get install mesa-common-dev  # optional, for OpenGL rendering
$ sudo apt-get install liblog4cxx10-dev
```
* Install the latest version of R (for Ubuntu 14.04) following the instruction at this [link](http://www.sysads.co.uk/2014/06/install-r-base-3-1-0-ubuntu-14-04/).
* Optional (but recommended): install the latest version of RStudio from this [link](http://www.rstudio.com/products/rstudio/download/).


### Mac OSX (Darwin)
* Install Xcode and XQuartz from the Apple site.
* Install the GSL libraries form source ([instructions](http://www.brianomeara.info/tutorials/brownie/gsl)). Alternatively it should be possible to install the binaries using Mac Ports.
* Install the log4cxx libraries from source [link](http://apache.fis.uniroma2.it/logging/log4cxx/0.10.0/apache-log4cxx-0.10.0.tar.gz ). (The installing procedure is similar to the GSL libraryies)
* Install the latest version of R from the official [web site](http://cran.rstudio.com/).
* Optional (but recommended): install the latest version of RStudio from this [link](http://www.rstudio.com/products/rstudio/download/).


### Rplanit installation
* To install the _Rplanit_ package from github (from the R prompt):
```
> install.packages("devtools") # it needs to be done only once
> library(devtools)
> install_github("planit-group/Rplanit", dependencies=TRUE)
```
* To load the package:
```
> library(Rplanit)
```


### Tutorial
A tutorial (in Italian) is available at http://rstudio-pubs-static.s3.amazonaws.com/66368_3aea5f348bd144ca89f9faf9bb1167b4.html
