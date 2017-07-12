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

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
# In this file functions to create artificial genotypes and phenotypes #
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#' Create a set of artificial genotypes
#'
#' This function return a matrix of 0/1/2... up to ploidy level.
#' Nothing fancy, each data point is sampled from an equally
#' probable distribution.
#'
#' @keywords internal
#'
#' @param samples number of individuals (rows)
#' @param markers number of SNPs (columns)
#' @param ploidy number of alleles per SNP
#'
#' @return a samples x markers matrix
artificialGenotypes = function(samples=50, markers=200, ploidy=2){
  genos = matrix(
    sample(x=ploidy + 1, size=samples * markers, replace = TRUE),
    nrow = samples, ncol = markers)
  return(genos - 1)
}

#' Create a set of artificial phenotypes
#'
#' Each marker is assigned an effect based on the passed distribution.
#' Heritability (h2) controls how much GEBV are close to phenotypes (h2=1 means
#' no distinction, h2=0 means no correlation)
#'
#' @param genotypes a samples x markers matrix of 0/1/2...ploidy
#' @param mu intercept (added to everything)
#' @param h2 heritability (must be greater than zero, an less than or equal to one)
#' @param markersEffectDistr name of the random generation function of
#'                           the selected distribution for markers effects. It must
#'                           accept n as argument indicating the number of elements
#'                           to be returned.
#' @param ... further parameters are passed to markersEffectDistr
#'
#' @keywords internal
#'
#' @return a list of three elements: $GEBV is an array of genetic breeding values,
#'        $phenotypes is the array of phenotypes, and $markerEffects is an array of
#'        marker effects
artificialPhenotypes = function(genotypes, mu=0, h2=0.8, markersEffectDistr = 'rnorm', ...){

  #marker effects
  parms = list(...)
  parms$n = ncol(genotypes)
  u = do.call(markersEffectDistr, parms)

  #computing GEBV and phenotypes
  GEBV = genotypes %*% u + mu
  phenos = GEBV + rnorm(nrow(genotypes),mean=0,sd=sqrt((1-h2)/h2*var(GEBV)))
  phenos = phenos[,1]

  #we are done
  return(list(
    GEBV = GEBV,
    phenotypes = phenos,
    markerEffects = u
  ))
}
