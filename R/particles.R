# PARTICLES --------------------------------------------------------------------
# set funzioni per la gestione/calcolo di quantità fisiche legate alle particelle


#' Z from ion
#'
#' Returns the charge of the atomic nucleus from the name of the ion. The ion can be specified from its element symbol, optionally including
#' its atomic mass number (e.g.: "C" or "12C" for carbon ion.)
#' Note: the available range is from 1H to 20 Ne.
#'
#' @param ions A character vector containing the ion names.
#' @param with.atomic.mass Uses ion names including the atomic mass number (default FALSE.)
#' @return A numeric vector containing the charges of the ions.
#'
#' @family Particles
#' @export
#'
z.from.ion <- function(ions, with.atomic.mass = FALSE) {                                                                                               
  if(with.atomic.mass) {ion.sequence <- c('1H', '4He', '7Li', '9Be', '11B', '12C', '14N', '16O', '19F', '20Ne')}                                       
  else {ion.sequence <- c('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne')}                                                                       
  match(ions, ion.sequence)                                                                                                                            
}

#' A from Ion
#'
#' Returns the atomic mass number (A) of the atomic nucleus from the name of the atom/ion.
#' The ion can be specified from its element symbol, optionally including
#' its atomic mass number (e.g.: "C" or "12C" for carbon ion.)
#' Note: the available range is from 1H to 20 Ne.
#'
#' @param ions A character vector containing the ion names.
#' @param with.atomic.mass Uses ion names including the atomic mass number (default FALSE.)
#' @return A numeric vector containing the atomic mass numbers.
#'
#' @family Particles
#' @export
#'
a.from.ion <- function(ions, with.atomic.mass = FALSE) {                                                                                               
  A.table <- c(1,4,7,9,11,12,14,16,19,20)                                                                                                              
  if(with.atomic.mass) {ion.sequence <- c('1H', '4He', '7Li', '9Be', '11B', '12C', '14N', '16O', '19F', '20Ne')}                                       
  else {ion.sequence <- c('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne')}                                                                       
  A.table[match(ions, ion.sequence)]                                                                                                                   
}


#' Ion from Z
#'
#' Returns the name of the ion/atom from its charge (Z).
#' The name is its element symbol, optionally including
#' its atomic mass number (e.g.: "C" or "12C" for carbon ion.)
#' Note: the available range is from 1H to 20 Ne.
#'
#' @param Z A numeric vector containing the ion charges.
#' @param with.atomic.mass Uses ion names including the atomic mass number (default FALSE.)
#' @return A numeric vector containing the atomic mass numbers.
#'
#' @family Particles
#' @export
#'
ion.from.z <- function(Z, with.atomic.mass = FALSE) {                                                                                                  
  Z[Z<0] <- NA                                                                                                                                         
  Z[Z>11] <- 11                                                                                                                                        
  if(with.atomic.mass) {ion.sequence <- c('1H', '4He', '7Li', '9Be', '11B', '12C', '14N', '16O', '19F', '20Ne', 'Z>10')}                               
  else {ion.sequence <- c('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Z>10')}                                                               
  ion.sequence[Z]                                                                                                                                      
}


#' Remove atomic mass number
#'
#' Returns the name of the ion/atom removing its atomic mass number specification (e.g. it returns "C" from "12C").
#' Note: the available range is from 1H to 20 Ne.
#'
#' @param ions A character vector containing the ion names (including the atomic mass number specification).
#' @return A character vector containing the ion names (without the atomic mass number specification).
#'
#' @family Particles
#' @export
#'
remove.atomic.mass <- function(ions) {                                                                                                                  
  ion.table <- c('1H', '4He', '7Li', '9Be', '11B', '12C', '14N', '16O', '19F', '20Ne')                                                                 
  ion.sequence <- c('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne')                                                                              
  ion.sequence[match(ions, ion.table)]                                                                                                                  
} 