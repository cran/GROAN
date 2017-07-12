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

#In this file we have functions and definitions necessary for workbench
#creation and functioning

####################### MAIN CLASS DEFINITION #######################
#' Workbench constructor
#'
#' This function creates a GROAN.Workbench instance (or fails trying). The created object contains:\cr
#' a) one regressor with its own specific configuration\cr
#' b) the crossvalidation parameters, describing the experiment.\cr
#' You can have a general description of the created object using the overridden \link{print.GROAN.Workbench}
#' function.\cr
#' It is possible to add other regressors to the created \code{GROAN.Workbench} object using \link{addRegressor}.
#' Once the \code{GROAN.Workbench} is created it must be passed to \link{GROAN.run} to start the experiment.\cr
#'
#' @param folds defaults to 10, as in 10-folds crossvalidation
#' @param reps number of times the whole test must be repeated, defaults to 5
#' @param stratified boolean indicating whether the crossvalidation should be standard (the default)
#'                   or stratified (keeping the same proportion of data strata in each fold). If no
#'                   strata are present in the \link[=createNoisyDataset]{GROAN.NoisyDataSet} object and
#'                   \code{stratified==TRUE} all samples will be considered belonging to the same strata ("dummyStrata")
#' @param outfolder folder where to save the data. If \code{NULL} (the default)
#'                  nothing will be saved. Filenames are standardized. If existing,
#'                  accuracy and hyperparameter files will be updated, otherwise are created. ExtraData
#'                  cannot be updated, so unique filenames will be generated using runId (see \link{GROAN.run})
#' @param saveHyperParms boolean indicating if the hyperparameters from regressor training should be
#'                       saved in \code{outfolder}. Defaults to FALSE.
#' @param saveExtraData boolean indicating if extradata from regressor training should be
#'                       saved in \code{outfolder} as R objects (using the \link[base]{save} function). Defaults to FALSE.
#' @param regressor regressor function. Defaults to \code{\link{phenoRegressor.rrBLUP}}
#' @param regressor.name string that will be used in reports. Keep that in mind when deciding names. Defaults to "default regressor"
#' @param ... extra parameter are passed to regressor function
#'
#' @return An instance of GROAN.Workbench
#' @seealso \link{addRegressor} \link{GROAN.run} \link{createNoisyDataset}
#' @examples #creating a Workbench with all default arguments
#' wb1 = createWorkbench()
#' #another Workbench, with different crossvalidation
#' wb2 = createWorkbench(folds=5, reps=20)
#' #a third one, with a different regressor and extra parameters passed to regressor function
#' wb3 = createWorkbench(regressor=phenoRegressor.BGLR, regressor.name='Bayesian Lasso', type='BL')
#' @export
createWorkbench = function (
  folds=10, reps=5, stratified=FALSE,
  outfolder=NULL, saveHyperParms=FALSE, saveExtraData=FALSE,
  regressor=phenoRegressor.rrBLUP,
  regressor.name = 'default regressor', ...){

  #checking inputs
  check.folds(folds)
  check.reps(reps)
  check.stratified(stratified)
  check.outfolder(outfolder)
  check.saves(saveHyperParms, saveExtraData, outfolder)
  check.regressor(regressor)
  check.regressor.name(regressor.name)

  #if we get here all parameters are good and we can go on

  #preparing the list of extra parms. Since there can be more
  #than one regressor, we create a list of list, using ordinals
  #as keys. So, given a GROAN.Workbench instance wb, all parameters
  #for first regressor are in wb$extraParms[[1]], for second
  #regressor are in wb$extraParms[[2]], and so forth. Actually,
  #because of list syntax, to obtain the original list a further
  #indexing is required, as in wb$extraParms[[1]][[1]], wb$extraParms[[2]][[1]],
  #and so forth.
  #Regressor functions are in wb$regressor[1], wb$regressor[2], and so on.
  #Regressor names (to be used in reports) aer in wb$regressor.name[1], wb$regressor.name[2] and so on.
  p = list()
  p[[1]] = list(...)

  #class attributes
  me = list(
    folds = folds,
    reps = reps,
    stratified = stratified,
    outfolder = outfolder,
    saveHyperParms = saveHyperParms,
    saveExtraData = saveExtraData,
    regressor = list(regressor),
    regressor.name = regressor.name,
    extraParms = p
  )

  ## Set the name for the class
  class(me) = c("GROAN.Workbench", class(me))

  #we are done
  return(me)
}

####################### INPUT CHECKERS #######################
#fails if folds is not a positive integer greater than 1
check.folds = function(folds){
  if (!is.naturalnumber(folds, 1)){
    stop(paste('Passed number of folds is not a positive integer greater than one:', folds), call. = FALSE)
  }
}

#fails if reps is not a positive integer
check.reps = function(reps){
  if (!is.naturalnumber(reps)){
    stop(paste('Passed number of repetitions is not a positive integer:', reps), call. = FALSE)
  }
}

#fails if outfolder is not valid string, or NULL
check.outfolder = function(outfolder){
  if (is.null(outfolder)){
    return()
  }

  #if we get here, is not NULL
  if (!is.string(outfolder)){
    stop(paste('Passed outfolder is not a valid string:', outfolder), call. = FALSE)
  }
}

#fails if either saveHyperParms or saveExtraData is not a boolean, or
#if it either is TRUE but no outfolder is specified
check.saves = function(saveHyperParms, saveExtraData, outfolder){
  #hyperparameters
  if (!is.boolean(saveHyperParms)){
    stop(paste('Passed saveHyperParms is not a valid boolean:', saveHyperParms), call. = FALSE)
  }
  if (saveHyperParms & is.null(outfolder)){
    stop(paste('Asked to save hyperparameters but no outfolder specified'), call. = FALSE)
  }

  #extra data
  if (!is.boolean(saveExtraData)){
    stop(paste('Passed saveExtraData is not a valid boolean:', saveExtraData), call. = FALSE)
  }
  if (saveExtraData & is.null(outfolder)){
    stop(paste('Asked to save saveExtraData but no outfolder specified'), call. = FALSE)
  }
}

#fails if stratified is not a valid boolean
check.stratified = function(stratified){
  if (!is.boolean(stratified)){
    stop(paste('Passed stratified is not a valid boolean:', stratified), call. = FALSE)
  }
}

#fails if regressor is not a valid function
check.regressor = function(regressor){
  if (!is.function(regressor)){
    stop(paste('Passed regressor is not a valid function. It\'s actually of type', typeof(regressor)), call. = FALSE)
  }

  #checking at least function arguments (there's no way to check return
  #values...)
  args.list = names(formals(regressor))
  if (!("phenotypes" %in% args.list) |
      !("genotypes" %in% args.list) |
      !("covariances" %in% args.list) |
      !("extraCovariates" %in% args.list)
      ){
    stop(paste('Passed regressor function must (at least) accept arguments  "phenotypes", "genotypes", "covariances" and "extraCovariates"'), call. = FALSE)
  }
}

#fails if regressor.name is not a string
check.regressor.name = function(regressor.name){
  if (!is.string(regressor.name)){
    stop(paste('Passed regressor.name is not a valid string:', regressor.name), call. = FALSE)
  }
}

#fails if the regressor name is already in the char array of regressor names
check.regressor.name.is.new = function(regressor.name, oldnames){
  if (regressor.name %in% oldnames){
    stop(paste('Duplicated regressor name:', regressor.name), call. = FALSE)
  }
}

####################### METHODS #######################
#' Add an extra regressor to a Workbench
#'
#' This function adds a regressor to an existing \link[=createWorkbench]{GROAN.Workbench} object.
#'
#' @param wb the GROAN.Workbench instance to be updated
#' @param regressor regressor function
#' @param regressor.name string that will be used in reports. Keep in mind that when deciding names.
#' @param ... extra parameters are passed to the regressor function
#'
#' @return an updated instance of the original GROAN.Workbench
#' @seealso \link{createWorkbench} \link{GROAN.run}
#' @examples #creating a Workbench with all default arguments
#' wb = createWorkbench()
#' #adding a second regressor
#' wb = addRegressor(wb, regressor = phenoRegressor.dummy, regressor.name = 'dummy')
#'
#' \dontrun{
#' #trying to add again a regressor with the same name would result in a naming conflict error
#' wb = addRegressor(wb, regressor = phenoRegressor.dummy, regressor.name = 'dummy')}
#' @export
addRegressor = function(wb, regressor, regressor.name=regressor, ...){
  #checking if it's a regressor
  check.regressor(regressor)
  check.regressor.name(regressor.name)

  #we also need to check that the regressor name is something new
  check.regressor.name.is.new(regressor.name, wb$regressor.name)

  #new regressor requires a new position
  p = length(wb$regressor) + 1

  #if we get here, we can add the passed regressor and parameters
  wb$regressor[[p]] = regressor
  wb$regressor.name = c(wb$regressor.name, regressor.name)
  wb$extraParms[[p]] = list(...)

  #and we are done
  return(wb)
}

#' Measure Performance of a Prediction
#'
#' This method returns several performance metrics for the passed
#' predictions.
#'
#' @param truevals true values
#' @param predvals predicted values
#'
#' @return A named array with the following fields:
#' \describe{
#'  \item{pearson}{Pearson's correlation}
#'  \item{spearman}{Spearmans' correlation (order based)}
#'  \item{rmse}{Root Mean Square Error}
#'  \item{mae}{Mean Absolute Error}
#'  \item{coeff_det}{Coefficient of determination}
#' }
#' @export
measurePredictionPerformance = function(truevals, predvals){
  #correlations
  pea = cor(truevals, predvals, method = 'pearson')
  spe = cor(truevals, predvals, method = 'spearman')

  #coeffiecient of determination
  ss_tot = sum((truevals - mean(truevals))^2)
  ss_res = sum((truevals - predvals)^2)
  r_square = 1 - ss_res / ss_tot

  return(c(
    pearson = pea,
    spearman = spe,
    cor_success = !is.na(pea),
    rmse = sqrt(mean((truevals - predvals)^2)),
    mae = mean(abs(truevals - predvals)),
    coeff_det = r_square
  ))
}

#' Print a GROAN Workbench object
#'
#' Short description for class GROAN.Workbench, created with \link{createWorkbench}.
#'
#' @param x object of class GROAN.Workbench.
#' @param ... ignored, put here to match S3 function signature
#'
#' @return This function returns the original \code{GROAN.Workbench} object invisibly (via \link[=invisible]{invisible(x)})
#' @export
print.GROAN.Workbench = function(x, ...){
  writeLines(paste(sep='', '--- EXPERIMENT ---'))
  writeLines(paste(sep='', 'folds          : ', x$folds))
  writeLines(paste(sep='', 'reps           : ', x$reps))
  writeLines(paste(sep='', 'stratified     : ', x$stratified))
  if (is.null(x$outfolder)){
    writeLines(paste(sep='', 'outfolder      : NULL'))
  }else{
    writeLines(paste(sep='', 'outfolder      : ', x$outfolder))
  }
  writeLines(paste(sep='', 'saveHyperParms : ', x$saveHyperParms))
  writeLines(paste(sep='', 'saveExtraData  : ', x$saveExtraData))

  writeLines(paste(sep='', '--- REGRESSORS ---'))
  for (i in 1:length(x$regressor.name)){
    writeLines(paste(sep='', 'Regressor ', i, '    : ', x$regressor.name[i]))
  }
  return(invisible(x))
}
