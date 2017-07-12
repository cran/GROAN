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

####################################################################################
# In this file we keep functions to add noise (of several types) on a phenotype    #
# array. The general form would be:                                                #
#                                                                                  #
# noiseInjector.distributionName = function(phenotypes, otherparams...)            #
#                                                                                  #
# where:                                                                           #
# - distributionName : name of the distribution from which the noise is sampled or #
#                      some other string describing the type of noise              #
# - phenotypes       : array of numbers                                            #
# - otherparams      : zero or more extra parameters necessary to qualify the      #
#                      noise injection operation(s)                                #
####################################################################################

#' Inject Normal Noise
#'
#' This function adds to the passed \code{phenotypes} array noise sampled from
#' a normal distribution with the specified mean and standard deviation.\cr
#' The function can interest the totality of the passed phenotype array or
#' a random subset of it (commanded by \code{subset} parameter).
#'
#' @param phenotypes an array of numbers.
#' @param mean mean of the normal distribution.
#' @param sd standard deviation of the normal distribution.
#' @param subset integer in [0,1], the proportion of original dataset to be injected
#'
#' @return An array, of the same size as phenotypes, where normal noise has been added to
#'         the original phenotype values.
#'
#' @examples
#' #a sinusoid signal
#' phenos = sin(seq(0,5, 0.1))
#' plot(phenos, type='p', pch=16, main='Original (black) vs. Injected (red), 100% affected')
#'
#' #adding normal noise to all samples
#' phenos.noise = noiseInjector.norm(phenos, sd = 0.2)
#' points(phenos.noise, type='p', col='red')
#'
#' #adding noise only to 30% of the samples
#' plot(phenos, type='p', pch=16, main='Original (black) vs. Injected (red), 30% affected')
#' phenos.noise.subset = noiseInjector.norm(phenos, sd = 0.2, subset = 0.3)
#' points(phenos.noise.subset, type='p', col='red')
#' @family noiseInjectors
#' @export
noiseInjector.norm = function(phenotypes, mean=0, sd=1, subset=1){
  #input check on subset, to avoid error
  if(!(is.single.slot(subset) & is.in.range(subset))){
    stop(paste('Passed subset is not in [0,1], received value:', subset), call. = FALSE)
  }

  #quick check
  if(subset == 0){
    #no noise added
    return(phenotypes)
  }

  #to number of involved samples
  num = ceiling(length(phenotypes) * subset)

  #selected at random
  sel = sample(length(phenotypes))[1:num]

  #creating the noise
  noise = rnorm(n=num,mean, sd)

  #injecting & returning
  phenotypes[sel] = phenotypes[sel] + noise
  return(phenotypes)
}

#' Inject Uniform Noise
#'
#' This function adds to the passed \code{phenotypes} array noise sampled from
#' a uniform distribution with the specified range.\cr
#' The function can interest the totality of the passed phenotype array or
#' a random subset of it (commanded by \code{subset} parameter).
#'
#' @param phenotypes an array of numbers.
#' @param min,max lower and upper limits of the distribution. Must be finite.
#' @param subset integer in [0,1], the proportion of original dataset to be injected
#' @return An array, of the same size as phenotypes, where uniform noise has been added to
#'         the original phenotype values.
#'
#' @examples
#' #a sinusoid signal
#' phenos = sin(seq(0,5, 0.1))
#' plot(phenos, type='p', pch = 16, main='Original (black) vs. Injected (red), 100% affected')
#'
#' #adding normal noise to all samples
#' phenos.noise = noiseInjector.unif(phenos, min=0.1, max=0.3)
#' points(phenos.noise, type='p', col='red')
#'
#' #adding noise only to 30% of the samples
#' plot(phenos, type='p', pch = 16, main='Original (black) vs. Injected (red), 30% affected')
#' phenos.noise.subset = noiseInjector.unif(phenos, min=0.1, max=0.3, subset = 0.3)
#' points(phenos.noise.subset, type='p', col='red')
#' @family noiseInjectors
#' @export
noiseInjector.unif = function(phenotypes, min=0, max=1, subset=1){
  #input check on subset, to avoid error
  if(!(is.single.slot(subset) & is.in.range(subset))){
    stop(paste('Passed subset is not in [0,1], received value:', subset), call. = FALSE)
  }

  #quick check
  if(subset == 0){
    #no noise added
    return(phenotypes)
  }

  #to number of involved samples
  num = ceiling(length(phenotypes) * subset)

  #selected at random
  sel = sample(length(phenotypes))[1:num]

  #creating the noise
  noise = runif(n=num, min, max)

  #injecting & returning
  phenotypes[sel] = phenotypes[sel] + noise
  return(phenotypes)
}

#' Noise Injector dummy function
#'
#' This noise injector does not add any noise. Passed \code{phenotypes} are
#' simply returned. This function is useful when comparing different
#' regressors on the same dataset without the effect of extra injected noise.
#'
#' @param phenotypes input phenotypes. This object will be returned without checks.
#'
#' @return the same passed \code{phenotypes}
#' @family noiseInjectors
#' @export
#' @examples
#' phenos = rnorm(10)
#' all(phenos == noiseInjector.dummy(phenos)) #TRUE
noiseInjector.dummy = function(phenotypes){
  return(phenotypes)
}

#' Swap phenotypes between samples
#'
#' This function introduces swap noise, i.e. a number of couples of
#' samples will have their phenotypes swapped.\cr
#' The number of couples is computed so that the total fraction of
#' interested phenotypes approximates \code{subset}.
#'
#' @param phenotypes an array of numbers
#' @param subset fraction of phenotypes to be interested by noise.
#'
#' @return the same passed \code{phenotypes}, but with some elements swapped
#' @examples
#' #a set of phenotypes
#' phenos = 1:10
#' #swapping two elements
#' phenos.sw2 = noiseInjector.swapper(phenos, 0.2)
#' #swapping four elements
#' phenos.sw4 = noiseInjector.swapper(phenos, 0.4)
#' #swapping four elements again, since 30% of 10 elements
#' #is rounded to 4 (two couples)
#' phenos.sw4.again = noiseInjector.swapper(phenos, 0.3)
#' @family noiseInjectors
#' @export
noiseInjector.swapper = function(phenotypes, subset = 0.1){
  #from subset to number of involved couples
  n = ceiling(length(phenotypes) * subset / 2)

  #if the number of couples is zero we return right away
  if (n==0){
    return(phenotypes)
  }

  #couples are created at random
  r = sample(length(phenotypes))

  #each couple is swapped
  for (i in 1:n){
    #chosing the samples to be swapped
    place1 = r[2*i-1]
    place2 = r[2*i]

    #swapping
    x = phenotypes[place1]
    phenotypes[place1] = phenotypes[place2]
    phenotypes[place2] = x
  }

  return(phenotypes)
}

