# Rplanit

RT Ion planning simulations and analysis

*Author*: Andrea Attili

*Maintainer*: Andrea Attili <attili@to.infn.it>

*Description*: A collection of methods and procedures to perform simulations and
analysis related to radiotherapy treatment planning (with ion beams). It includes a set
of radiobiological evaluations and display possibilities.


## Installation

### Ubuntu (12.04 64bit)
* The following libraries need to be installed from the Ubuntu repository:
```
$ sudo apt-get install libcurl4-openssl-dev
$ sudo apt-get install tk8.5-dev
$ sudo apt-get install libsdl1.2-dev
$ sudo apt-get install mesa-common-dev
$ sudo apt-get install libssl-dev
$ sudo ln -s /lib/x86_64-linux-gnu/libcrypto.so.1.0.0 /lib/x86_64-linux-gnu/libcrypto.so.6
$ sudo ln -s /lib/x86_64-linux-gnu/libssl.so.1.0.0 /lib/x86_64-linux-gnu/libssl.so.6
$ sudo apt-get install libgsl0-dev
$ sudo apt-get install liblog4cxx10-dev
```
* Install the latest version of R (for Ubuntu 12.04) following the instruction at this [link](http://livesoncoffee.wordpress.com/2012/12/09/installing-r-on-ubuntu-12-04/).
* Optional (but recommended): install the latest version of RStudio from this [link](http://www.rstudio.com/products/rstudio/download/).


### Mac OSX (Darwin)
* Install Xcode and XQuartz from the Apple site.
* Install the GSL libraries form source ([instructions](http://www.brianomeara.info/tutorials/brownie/gsl)). Alternatively it should be possible to install the binaries using Mac Ports.
* Install the latest version of R from the official [web site](http://cran.rstudio.com/).
* Optional (but recommended): install the latest version of RStudio from this [link](http://www.rstudio.com/products/rstudio/download/).


### Rplanit installation
* To install the _Rplanit_ package from github (from the R prompt):
```
> install.packages("devtools") # it needs to be done only once
> library(devtools)
> install_github("Rplanit", username="planit-group", dependencies=TRUE)
```
* To load the package:
```
> library(Rplanit)
```
* To install _PlanKIT_ (the Treatment Planning System kernel), _Gate_ (to perform Monte Carlo simulations) and _Survival_ (the code for the radiobiological simulations):
```
> install.plankit(<lincence number>) # please read the note below.
> install.gate()
> instal.surival()
> # et voila
```
* To update everything to the latest version, the last commands should be runned again.

#### Important Note
While _Rplanit_, _Gate_ and _Survivals_ are open-source, the TPS kernel, _PlanKIT_ is proprietary by INFN/IBA. To use it for research purposes please contact Andrea Attili (attili@to.infn.it) to obtain the licence number and have access to the binaries.