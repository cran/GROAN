# Copyright 2017 Nelson Nazzicari
# This file is part of GROAN.
#
# GROAN is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GROAN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GROAN  If not, see <http://www.gnu.org/licenses/>.

##########################################################
# In this script, common functionalities and handy hacks #
##########################################################

#' Is a positive integer?
#'
#' \code{is.naturalnumber} returns true if the passed argument
#' is a positive integer, false otherwise.
#' Implementation taken from Marcog's answer
#' \href{http://stackoverflow.com/questions/4562257/what-is-the-fastest-way-to-check-if-a-number-is-a-positive-natural-number-in-r}{to this question}
#'
#' @param x the value(s) to be tested
#' @param low.Limit the greatest value not accepted. Defaults to zero, meaning that
#'                 one is the smallest integer that returns true.
#'
#' @keywords internal
#'
#' @return TRUE if the value is a positive integer, FALSE otherwise (or NA)
is.naturalnumber = function(x, low.Limit=0) {
  #base validity check
  if (any(is.na(x))) return (FALSE)
  if (any(is.null(x))) return (FALSE)

  #numeric check
  tol = .Machine$double.eps^0.5
  return (x > low.Limit & abs(x - round(x)) < tol)
}

#' Is a string?
#'
#' Returns TRUE if the passed x variable is a length one variable
#' containing characters values
#'
#' @param x the variable to be checked
#'
#' @keywords internal
#'
#' @return a boolean
is.string = function(x){
  return (length(x) == 1 & all(is.character(x)))
}


#' Is a boolean?
#'
#' Returns TRUE if the passed x variable is a length one variable
#' containing a valid TRUE/FALSE value.
#'
#' @param x the variable to be checked
#'
#' @keywords internal
#'
#' @return a boolean
is.boolean = function(x){
  if (length(x) > 1){
    return (FALSE)
  }

  return (x == TRUE | x == FALSE)
}

#create a temporary directory. Arguments are passed to tempfile()
tempdir.create = function(...){
  dir = tempfile(...)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  return(dir)
}

#returns a log header, useful for separating log entries
#or stuff that is appended. Passed text is put in the
#log header, too.
log.header = function(text = NULL){
  return(paste(sep="\n",
    "===================================",
    date(),
    text,
    "-----------------------------------",
    '' #this is needed to end with a newline
  ))
}

#write the passed list to the passed connection
#crashes on list of lists
writelist = function(x, connection, level=0){
  pr = paste(rep_len('-', level), collapse = '', sep='')
  for (i in 1:length(x)){
    #has this slot of the list a name?
    n = names(x[i])
    writeLines(paste(pr, '>[', i, '] ', n, collapse = '', sep=''), con = connection)
    if (is.list(x[[i]])){
      writelist(x[[i]], connection, level + 1)
    }else{
      if (typeof(x[[i]]) == 'language'){
        #objects of type "language" (typically the "call" fields)
        #require a special handling
        write(deparse(x[[i]], width.cutoff = 500), file = connection)
      }else{
        #standard object pose no problem. We just try to adjust the column width,
        #if it makes any sense
        ncolumns.here = ncol(x[[i]])

        #asking ncol to unknown object can lead to various results: NaN, NA, NA_real_, NULL...
        #we try to be conservative and cover all bases
        if (is.null(ncolumns.here)) ncolumns.here = 20
        if (!is.integer(ncolumns.here) | is.na(ncolumns.here))  ncolumns.here = 20

        #we can finally write
        write(x[[i]], file = connection, ncolumns = ncolumns.here)
      }
    }
  }
}

#' Assignments for stratified crossvalidation
#'
#' This function creates the assignments to perform a
#' stratified crossvalidation. It requires, for each
#' element of the set, a label describing the strata it
#' belongs to. It also requires the number of folds.
#'
#' A warning is triggered if the number of folds is greater than the number of
#' elements of any stratum.
#'
#' @keywords internal
#'
#' @param strata an array of strata (will be treated as a factor)
#' @param folds.num the number of folds
#'
#' @return an array, same length as \code{strata}, of numbers in the 1:\code{folds.num} set
stratified_crossvalidation_folds = function(strata, folds.num){
  res = rep(NA, length(strata))
  #for each different stratum
  for (c in unique(strata)){
    #where are placed the elements of this stratum
    where = strata == c

    #how many elements in this stratum
    n = sum(where)

    #small check: do we have enough elements for the folds
    if (n < folds.num){
      warning(paste('Stratum', c, 'has less elements than requested folds'))
    }

    #the folds for this stratum
    str = rep(1:folds.num, length.out = n)

    #shuffling
    str = str[sample(n)]

    #placing them in the original result
    res[where] = str
  }

  return(res)
}


#' Is in the passed numeric range
#'
#' Function returns TRUE if the passed variable is numeric and all its content
#' is in the passed range (defined by \code{min} and \code{max}). Works also in
#' array (in that case all values must be in the range).
#'
#' @param x value to be tested
#' @param min minimum admitted value
#' @param max maximum admitted value
#'
#' @return TRUE if x is a numeric completely in the specified range, FALSE otherwise
#' @keywords internal
#' @export
is.in.range = function(x, min=0, max=1){
  if (is.numeric(x)){
    if (all(x <= max)){
      if(all(x >= min)){
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

#' Is a single slot thing?
#'
#' Checks if the passed variable is a single slot thing, meaning it contains only a single
#' value (numeric, character, whatever) and nothing more.
#' Works with array, vectors, matrix, data.frame...
#'
#' @param x the thing to be tested
#' @param NULL.is.single should NULL be considered a single slot or not (default: not)
#'
#' @return TRUE if is single slot, FALSE otherwise
#' @keywords internal
#' @export
#'
#' @examples
#' is.single.slot(5)   #TRUE
#' is.single.slot('foobar')   #TRUE
#' is.single.slot(NULL)       #depends on NULL.is.single
#' is.single.slot(NA)         #TRUE
#' is.single.slot(c(1,2,5))   #FALSE
#' is.single.slot(matrix(0, 10, 5))   #FALSE
#' is.single.slot(matrix(0, 1, 1))   #TRUE
is.single.slot = function(x, NULL.is.single = FALSE){
  #NULL is a very special thing
  if (NULL.is.single){
    if (is.null(x)){
      return(TRUE)
    }
  }

  #checking for single dimensional objects (without the dim thing)
  if (is.null(dim(x))){
    return(length(x) == 1)
  }

  #if we get here is multidimensional, let's find its max span
  return(max(dim(x)) == 1)
}
