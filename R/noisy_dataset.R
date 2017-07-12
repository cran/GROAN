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

#in this file functions related to creation of noisy dataset

####################### MAIN CLASS DEFINITION #######################
#' Noisy Data Set Constructor
#'
#' This function creates a GROAN.NoisyDataset object (or fails trying). The
#' class will contain all noisy data set components: genotypes and/or covariance matrix,
#' phenotypes, strata (optional), a noise injector function and its parameters.\cr
#' You can have a general description of the created object using the overridden \link{print.GROAN.NoisyDataset}
#' function.
#'
#' @param name A string defining the dataset name, used later do identify this
#'             particular instance in reports and result files. It is advisable for
#'             it to be it somewhat meaningful (to you, GROAN simply reports it as it is)
#' @param genotypes Matrix or dataframe containing SNP genotypes, one row per sample (N), one column per marker (M), 0/1/2 format (for diploids)
#'                  or 0/1/2.../ploidy in case of polyploids
#' @param covariance square (NxN) matrix of covariances between samples
#' @param phenotypes numeric array, M slots
#' @param strata array of M slots, describing the strata each data point belongs to. This is
#'                used for stratified crossvalidation (see \code{\link{createWorkbench}})
#' @param extraCovariates dataframe of optional extra covariates (N lines, one column per extra covariate).
#'                        Numeric ones will be normalized, string and categorical ones will be transformed
#'                        in stub TRUE/FALSE variables (one per possible value, see \link[stats]{model.matrix}).
#' @param ploidy number of haploid sets in the cell. Defaults to 2 (diploid).
#' @param noiseInjector name of a noise injector function, defaults to \link{noiseInjector.dummy}
#' @param ... further arguments are passed along to noiseInjector
#'
#' @return a GROAN.NoisyDataset object.
#' @seealso \link{GROAN.run} \link{createNoisyDataset}
#' @export
#'
#' @examples #For more complete examples see the package vignette
#' #creating a noisy dataset with normal noise
#' nds = createNoisyDataset(
#'   name = 'PEA, normal noise',
#'   genotypes = GROAN.pea.SNPs,
#'   phenotypes = GROAN.pea.yield,
#'   noiseInjector = noiseInjector.norm,
#'   mean = 0,
#'   sd = sd(GROAN.pea.yield) * 0.5
#' )
createNoisyDataset = function (name, genotypes=NULL, covariance=NULL, phenotypes, strata=NULL, extraCovariates=NULL, ploidy=2, noiseInjector=noiseInjector.dummy, ...){
  #checking inputs
  check.name(name)
  check.ploidy(ploidy)
  check.phenotypes(phenotypes)
  check.genetic_data(genotypes, covariance, n=length(phenotypes), ploidy=ploidy)
  check.strata(strata, n=length(phenotypes))
  check.extraCovariates(extraCovariates, n=length(phenotypes))
  check.noiseInjector(noiseInjector)

  #class attributes
  me = list(
    name = name,
    ploidy = ploidy,
    genos = genotypes,
    covs = covariance,
    phenos = phenotypes,
    strata = strata,
    injector = noiseInjector,
    extraParms = list(...)
  )

  #optional extra covariates get special treatment
  if (is.null(extraCovariates)){
    me$extraCovs = NULL
  }else{
    #if extraCovs is already a multi dimensional thing (matrix, dataframe, etc)
    #we are good (correct dimentions, we have colnames). Otherwise
    #we force it
    if (is.null(dim(extraCovariates))){
      extraCovariates = as.data.frame(extraCovariates)
    }
    if (length(dim(extraCovariates)) == 1){
      extraCovariates = as.data.frame(extraCovariates)
    }

    me$extraCovs = model.matrix( ~ ., extraCovariates)
    #because of how we built the model matrix, the first
    #column is "(Intercept)". This is ok, but the parenthesis
    #in the column name are bad for some regressors. Let's
    #correct it
    tmp.names = colnames(me$extraCovs)
    tmp.names[1] = 'Intercept'
    colnames(me$extraCovs) = tmp.names

    #Let's just keep the original covariate names
    me$extraCovsOriginalNames = colnames(extraCovariates)
  }

  ## Set the name for the class
  class(me) = c("GROAN.NoisyDataset", class(me))

  #we are done
  return(me)
}

####################### INPUT CHECKERS #######################
#fails if phenotypes is not numeric, or if it contains NA
check.phenotypes = function(phenotypes){
  if(!is.numeric(phenotypes)){
    stop(paste('Passed phenotypes is not a numeric array.'), call. = FALSE)
  }
  if(any(is.na(phenotypes))){
    stop(paste('Passed phenotypes contains missing values.'), call. = FALSE)
  }
}

#fails if ploidy is not a positive integer
check.ploidy = function(ploidy){
  if (!is.naturalnumber(ploidy)){
    stop(paste('Passed ploidy is not a positive integer:', ploidy), call. = FALSE)
  }
}

#check that extra covariates are NULL or a matrix/df in the form:
#-n rows (where n is the number of samples)
#-no missings
#fails if any condition is not met
check.extraCovariates = function(extraCovariates, n){
  #extraCovariates are optional, so NULL is a valid value
  if (is.null(extraCovariates)) return()

  #adjusting for single dimension (vector) vs. multiple dimension (array, data.frame, matrix)
  if (is.null(dim(extraCovariates))){
    extraCovariates = as.data.frame(extraCovariates)
  }
  if (length(dim(extraCovariates)) == 1){
    extraCovariates = as.data.frame(extraCovariates)
  }

  #dimensional check
  if (nrow(extraCovariates) != n){
    stop(paste('Passed extra covariates should have as many rows as phenotypes slots.'), call. = FALSE)
  }

  #missing values
  if (any(is.na(extraCovariates))){
    stop(paste('Passed extra covariates contains missing values.'), call. = FALSE)
  }
}

#checks that at least one of genotypes and covariance is not NULL
#calls the proper further checks on the available data
check.genetic_data = function(genotypes, covariance, n, ploidy){
  if (all(is.null(genotypes), is.null(covariance))){
    stop(paste('At least one among genotypes and covariance matrix need to be not NULL'), call. = FALSE)
  }
  if (!is.null(genotypes)){
    check.genotypes(genotypes, n, ploidy)
  }
  if (!is.null(covariance)){
    check.covariance(covariance, n)
  }
}

#check genotypes are a matrix/df in the form
#-n rows (where n is the number of samples)
#-one column per locus
#-values in 0/1/2/...ploidy
#-no missings
#fails if any condition is not met
check.genotypes = function(genotypes, n, ploidy){
  #dimensional check
  if (nrow(genotypes) != n){
    stop(paste('Passed genotypes should have as many rows as phenotypes slots.'), call. = FALSE)
  }

  #missing values
  if (any(is.na(genotypes))){
    stop(paste('Passed genotypes contains missing values.'), call. = FALSE)
  }

  #for numeric check we need a matrix
  genotypes.max = as.matrix(genotypes)

  #The following checks were removed to allow experimenting
  #with kinship values as genotypes. Please reactivate it
  #once kinship has proper support

  #are all valid integer values?
    if(!all(is.naturalnumber(genotypes.max, -1))){
      stop(paste('Passed genotypes contains non integer or negative values.'), call. = FALSE)
    }

  #are all in the ploidy range?
  if(max(genotypes.max > ploidy)){
   stop(paste('Passed genotypes contains value(s) greater than ploidy.'), call. = FALSE)
  }
}

#check genotypes are a matrix/df in the form
#-n x n (where n is the number of samples)
#-no missings
#fails if any condition is not met
check.covariance = function(covariance, n){
  #dimensional check
  if (nrow(covariance) != n){
    stop(paste('Passed covariance should have as many rows as phenotypes slots.'), call. = FALSE)
  }
  if (ncol(covariance) != n){
    stop(paste('Passed covariance should have as many columns as phenotypes slots.'), call. = FALSE)
  }

  #missing values
  if (any(is.na(covariance))){
    stop(paste('Passed covariance contains missing values.'), call. = FALSE)
  }
}

#check that strata array is either NULL or contain n elements,
#with no empty ones
check.strata = function(strata, n){
  #if it's NULL it's ok
  if(is.null(strata)){
    return()
  }

  #dimensional check
  if (length(strata) != n){
    stop(paste('Passed strata should have as many rows as phenotypes slots.'), call. = FALSE)
  }

  #no empty spaces allowed
  if(any(is.na(strata))){
    stop(paste('Passed strata contains missing values.'), call. = FALSE)
  }
}

#fails if noiseInjector is not a valid noise injection function or function name
check.noiseInjector = function(noiseInjector){
  if (!is.function(noiseInjector)){
    stop(paste('Passed noiseInjector is not a function'), call. = FALSE)
  }

  #checking at least function arguments (there's no way to check return
  #values...)
  args.list = names(formals(noiseInjector))
  if (!("phenotypes" %in% args.list)){
    stop(paste('Passed noiseInjector function must (at least) accept arguments "phenotypes"'), call. = FALSE)
  }
}

#fails if name is not a string
check.name = function(name){
  if (!is.string(name)){
    stop(paste('Passed name is not a valid string:', name), call. = FALSE)
  }
}

###################### METHODS #######################
#' Generate an instance of noisy phenotypes
#'
#' Given a \code{Noisy Dataset} object, this function
#' applies the noise injector to the data and returns
#' a noisy version of it.
#' It is useful for inspecting the noisy injector effects.
#' @param nds a \code{Noisy Dataset} object
#' @return the phenotypes contained in \code{nds} with added noise.
#'
#' @export
getNoisyPhenotype = function(nds){
  #the noise injection function
  f = nds$injector
  #the phenotypes
  p = nds$phenos
  #the extra arguments
  ex = nds$extraParms
  #creating an argument list where phenotypes are the first element
  parms = list(phenotypes=p)
  parms = c(parms, ex)

  #done
  return (do.call(f, parms))
}

#' String of extra covariates names
#'
#' Given a \code{Noisy Dataset} object, this function
#' returns a representation of the extra covariates
#' present in the object. If no extra covariates are
#' present an empty string is returned.
#'
#' @keywords internal
#'
#' @param nds a \code{Noisy Dataset} object.
#' @param separator used for string concatenation.
#' @return a string representation of extra covariates names.
getExtraCovariatesNames = function(nds, separator=' '){
  if (is.null(nds$extraCovs)){
    return('')
  }
  return(paste(collapse=separator, nds$extraCovsOriginalNames))
}

#' Print a GROAN Noisy Dataset object
#'
#' Short description for class GROAN.NoisyDataset, created with \link{createNoisyDataset}.
#'
#' @param x object of class GROAN.NoisyDataset.
#' @param ... ignored, put here to match S3 function signature
#'
#' @return This function returns the original \code{GROAN.NoisyDataset} object invisibly (via \link[=invisible]{invisible(x)})
#' @export
print.GROAN.NoisyDataset = function(x, ...){
  writeLines(paste(sep='', 'name              : ', x$name))
  writeLines(paste(sep='', 'ploidy            : ', x$ploidy))
  writeLines(paste(sep='', 'samples           : ', length(x$phenos)))

  SNPs.num = 'absent'
  if (!is.null(x$genos)){
    SNPs.num = ncol(x$genos)
  }
  writeLines(paste(sep='', 'SNPs              : ', SNPs.num))

  covariance = 'present'
  if (is.null(x$covs)){
    covariance = 'absent'
  }
  writeLines(paste(sep='', 'covariance matrix : ', covariance))

  extraCovariates = 'absent'
  if (!is.null(x$extraCovariates)){
    extraCovariates = paste(collapse = ', ', colnames(x$extraCovariates))
  }
  writeLines(paste(sep='', 'extra covariates  : ', extraCovariates))

  strata = 'absent'
  if (!is.null(x$strata)){
    strata = paste(collapse = ', ', unique(x$strata))
  }
  writeLines(paste(sep='', 'strata           : ', strata))

  return(invisible(x))
}
