# Script di battaglia per caricare le librerie al volo (per modifiche durante altre sessioni di lavoro)

# librerie esterne  
library(ggplot2)
library(fields)
library(splancs)
library(markdown)
library(rgl)
library(misc3d)
library(oro.nifti)
library(oro.dicom)
library(gtable)
library(knitr)
  
# carica tutte le funzioni definite nel folder...
  
# variabili di environment
my.folder  <- '/home/andrea/Programmi/Rplanit/R/'
 
my.files <- dir(my.folder, pattern=glob2rx('*.R'))
my.files <- my.files[my.files!='test-functions.R']
  
for(my.file in my.files) {
  source(paste(my.folder, my.file, sep=''))
}

