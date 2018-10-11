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

########################################################################################
# In this file functions to do the actual regression on phenotypes after a proper      #
# training. The functions would probably be wrappers around existing                   #
# functions/libraries more than actualy implementation of mathematical models.         #
# The general form is:                                                                 #
#                                                                                      #
# phenoRegressor.algo = function(                                                      #
#           phenotypes, genotypes, extraCovariates, covariances, otherparams)          #
#                                                                                      #
# where:                                                                               #
# - algo            : name of the regressor                                            #
# - phenotypes      : phenotypes, a numeric array (n x 1), missing values are          #
#                     to be predicted                                                  #
# - genotypes       : SNP genotypes, one row per phenotype (n), one column per         #
#                     marker (m), values in 0/1/2 for diploids or 0/1/2/...ploidy for  #
#                     polyploids. Can be NULL if \code{covariances} is present.        #
# - covariances     : square matrix (n x n) of covariances. Can be NULL if genotypes   #
#                     is present.                                                      #
# - extraCovariates : extra covariates set, one row per phenotype (n), one column per  #
#                     covariate (w). If NULL no extra covariates are considered.       #
# - otherparams     : zero or more extra parameters necessary to specify the behaviour #
#                     of the regression algorithm                                      #
#                                                                                      #
# Checks on input (type, size) are not required, since these functions are supposed to #
# be used in a GROAN experiment, where checks are done when NoisyDataset and Workbench #
# are created.                                                                         #
#                                                                                      #
# The function should return a list with the following fields:                         #
# - 'predictions' : an array of (n x 1) of phenotypes with all NAs filled (predicted)  #
# - 'hyperparams' : (optional) named array of hyperparameters selected during          #
#                   training                                                           #
# - 'extradata'   : (optional) whatever extra data/objects are produced during         #
#                   training                                                           #
########################################################################################

#' Regression dummy function
#'
#' This function is for development purposes. It returns, as "predictions", an array of
#' random numbers. It accept the standard inputs and produces a formally correct output. It
#' is, obviously, quite fast.
#'
#' @param phenotypes phenotypes, numeric array (n x 1), missing values are predicted
#' @param genotypes SNP genotypes, one row per phenotype (n), one column per marker (m), values in 0/1/2 for
#'                  diploids or 0/1/2/...ploidy for polyploids. Can be NULL if \code{covariances} is present.
#' @param covariances square matrix (n x n) of covariances. Can be NULL if \code{genotypes} is present.
#' @param extraCovariates extra covariates set, one row per phenotype (n), one column per covariate (w).
#'                 If NULL no extra covariates are considered.
#' @return The function should return a list with the following fields:
#' \itemize{
#'   \item \code{predictions} : an array of (k) predicted phenotypes
#'   \item \code{hyperparams} : named array of hyperparameters selected during training
#'   \item \code{extradata}   : any extra information
#' }
#'
#' @family phenoRegressors
#'
#' @examples #genotypes are not really investigated. Only
#' #number of test phenotypes is used.
#' phenoRegressor.dummy(
#'   phenotypes = c(1:10, NA, NA, NA),
#'   genotypes = matrix(nrow = 13, ncol=30)
#' )
#' @export
phenoRegressor.dummy = function (phenotypes, genotypes, covariances, extraCovariates){
  res = list(
    predictions = runif(length(phenotypes)),
    hyperparams = c(some = 1, params='a'),
    extradata = 'filler extra data for phenoRegressor.dummy'
  )
  return(res)
}

#' Regression using BGLR package
#'
#' This is a wrapper around \code{\link[BGLR]{BGLR}}. As such, it won't work if BGLR package
#' is not installed.\cr
#' Genotypes are modeled using the specified \code{type}. If \code{type} is 'RKHS' (and only
#' in this case) the covariance/kinship matrix \code{covariances} is required, and it will be modeled
#' as matrix K in BGLR terms. In all other cases genotypes and covariances are put in the model
#' as X matrices.\cr
#' Extra covariates, if present, are modeled as FIXED effects.
#'
#' @param phenotypes phenotypes, a numeric array (n x 1), missing values are predicted
#' @param genotypes SNP genotypes, one row per phenotype (n), one column per marker (m), values in 0/1/2 for
#'                  diploids or 0/1/2/...ploidy for polyploids. Can be NULL if \code{covariances} is present.
#' @param covariances square matrix (n x n) of covariances. Can be NULL if \code{genotypes} is present.
#' @param extraCovariates extra covariates set, one row per phenotype (n), one column per covariate (w).
#'                 If NULL no extra covariates are considered.
#' @param type character literal, one of the following: FIXED (Flat prior), BRR (Gaussian prior),
#'             BL (Double-Exponential prior), BayesA (scaled-t prior),
#'             BayesB (two component mixture prior with a point of mass at zero and a scaled-t slab),
#'             BayesC (two component mixture prior with a point of mass at zero and a Gaussian slab)
#' @param ... extra parameters are passed to \code{\link[BGLR]{BGLR}}
#'
#' @return The function returns a list with the following fields:
#' \itemize{
#'   \item \code{predictions} : an array of (n) predicted phenotypes, with NAs filled and all other positions repredicted (useful for calculating residuals)
#'   \item \code{hyperparams} : empty, returned for compatibility
#'   \item \code{extradata}   : list with information on trained model, coming from \code{\link[BGLR]{BGLR}}
#' }
#' @family phenoRegressors
#' @seealso \link[BGLR]{BGLR}
#' @export
#' @examples
#' \dontrun{
#' #using the GROAN.KI dataset, we regress on the dataset and predict the first ten phenotypes
#' phenos = GROAN.KI$yield
#' phenos[1:10]  = NA
#'
#' #calling the regressor with Bayesian Lasso
#' results = phenoRegressor.BGLR(
#'   phenotypes = phenos,
#'   genotypes = GROAN.KI$SNPs,
#'   covariances = NULL,
#'   extraCovariates = NULL,
#'   type = 'BL', nIter = 2000 #BGLR-specific parameters
#' )
#'
#' #examining the predictions
#' plot(GROAN.KI$yield, results$predictions,
#'      main = 'Train set (black) and test set (red) regressions',
#'      xlab = 'Original phenotypes', ylab = 'Predicted phenotypes')
#' points(GROAN.KI$yield[1:10], results$predictions[1:10], pch=16, col='red')
#'
#' #printing correlations
#' test.set.correlation  = cor(GROAN.KI$yield[1:10], results$predictions[1:10])
#' train.set.correlation = cor(GROAN.KI$yield[-(1:10)], results$predictions[-(1:10)])
#' writeLines(paste(
#'   'test-set correlation :', test.set.correlation,
#'   '\ntrain-set correlation:', train.set.correlation
#' ))
#' }
phenoRegressor.BGLR = function (phenotypes, genotypes, covariances, extraCovariates,
                                type = c('FIXED', 'BRR', 'BL', 'BayesA', 'BayesB', 'BayesC', 'RKHS'),
                                ...){
  #is BGLR installed?
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("BGLR package needed for this regressor to work. Please install it.",
         call. = FALSE)
  }

  #checking regression type in a supported list
  type = match.arg(type)

  #preparing the list to call the BGLR
  BGLR.args = list(...)
  BGLR.args$y = phenotypes

  #what type of ETA depends on the user's choice of model
  #in particular, RKHS needs a kinship/covariances, while all other can work on
  #either SNPs or covariances
  if (type == 'RKHS'){
    #this will fail if covariances are not present
    BGLR.args$ETA = list(list(K=covariances, model=type))
  }else{
    #if we have snps, we add the snps
    BGLR.args$ETA = list()
    if (!is.null(genotypes)){
      BGLR.args$ETA = list(list(X=genotypes, model=type))
    }
    #if we have kinship, we add kinship
    if (!is.null(covariances)){
      num = length(BGLR.args$ETA)
      BGLR.args$ETA[[num+1]] = list(X=covariances, model=type)
    }
  }

  #if we have extra covariates, we put them in as fixed effects
  if (!is.null(extraCovariates)){
    num = length(BGLR.args$ETA)
    BGLR.args$ETA[[num+1]] = list(X=extraCovariates, model='FIXED')
  }

  #since there is no way to NOT save the progress files
  #we create a tmpdir and then delete it at the end
  tmpdir = tempdir.create('BGLR')
  BGLR.args$saveAt = tmpdir

  #we should run with verbose=FALSE, but the user takes
  #precedence
  if (is.null(BGLR.args$verbose)){
    BGLR.args$verbose = FALSE
  }

  #ready to run
  trained.model = do.call(BGLR::BGLR, args = BGLR.args)

  #delete tempdir
  unlink(tmpdir, recursive = TRUE)

  #creating the result list
  return(list(
    predictions = trained.model$yHat,
    hyperparams = c(Note = 'No hyperparameters saved from BGLR, please see extra data'),
    extradata = trained.model
  ))
}

#' SNP-BLUP or G-BLUP using rrBLUP package
#'
#' This is a wrapper around \code{rrBLUP} function \code{\link[rrBLUP]{mixed.solve}}.
#' It can either work with genotypes (in form of a SNP matrix) or with kinships (in form of a covariance
#' matrix). In the first case the function will implement a SNP-BLUP, in the second a G-BLUP. An error is
#' returned if both SNPs and covariance matrix are passed.\cr
#' In rrBLUP terms, genotypes are modeled as random effects (matrix Z), covariances as matrix K, and
#' extra covariates, if present, as fixed effects (matrix X).\cr
#' Please note that this function won't work if rrBLUP package is not installed.
#'
#' @param phenotypes phenotypes, a numeric array (n x 1), missing values are predicted
#' @param genotypes SNP genotypes, one row per phenotype (n), one column per marker (m), values in 0/1/2 for
#'                  diploids or 0/1/2/...ploidy for polyploids. Can be NULL if \code{covariances} is present.
#' @param covariances square matrix (n x n) of covariances.
#' @param extraCovariates optional extra covariates set, one row per phenotype (n), one column per covariate (w).
#'                 If NULL no extra covariates are considered.
#' @param ... extra parameters are passed to rrBLUP::mixed.solve
#'
#' @return The function returns a list with the following fields:
#' \itemize{
#'   \item \code{predictions} : an array of (k) predicted phenotypes
#'   \item \code{hyperparams} : named vector with the following keys: Vu, Ve, beta, LL
#'   \item \code{extradata}   : list with information on trained model, coming from \code{\link[rrBLUP]{mixed.solve}}
#' }
#' @family phenoRegressors
#' @seealso \link[rrBLUP]{mixed.solve}
#' @export
#' @examples
#' \dontrun{
#' #using the GROAN.KI dataset, we regress on the dataset and predict the first ten phenotypes
#' phenos = GROAN.KI$yield
#' phenos[1:10]  = NA
#'
#' #calling the regressor with ridge regression BLUP on SNPs and kinship
#' results.SNP.BLUP = phenoRegressor.rrBLUP(
#'   phenotypes = phenos,
#'   genotypes = GROAN.KI$SNPs,
#'   SE = TRUE, return.Hinv = TRUE #rrBLUP-specific parameters
#' )
#' results.G.BLUP = phenoRegressor.rrBLUP(
#'   phenotypes = phenos,
#'   covariances = GROAN.KI$kinship,
#'   SE = TRUE, return.Hinv = TRUE #rrBLUP-specific parameters
#' )
#'
#' #examining the predictions
#' plot(GROAN.KI$yield, results.SNP.BLUP$predictions,
#'      main = '[SNP-BLUP] Train set (black) and test set (red) regressions',
#'      xlab = 'Original phenotypes', ylab = 'Predicted phenotypes')
#' abline(a=0, b=1)
#' points(GROAN.KI$yield[1:10], results.SNP.BLUP$predictions[1:10], pch=16, col='red')
#'
#' plot(GROAN.KI$yield, results.G.BLUP$predictions,
#'      main = '[G-BLUP] Train set (black) and test set (red) regressions',
#'      xlab = 'Original phenotypes', ylab = 'Predicted phenotypes')
#' abline(a=0, b=1)
#' points(GROAN.KI$yield[1:10], results.G.BLUP$predictions[1:10], pch=16, col='red')
#'
#' #printing correlations
#' correlations = data.frame(
#'   model = 'SNP-BLUP',
#'   test_set_correlations = cor(GROAN.KI$yield[1:10], results.SNP.BLUP$predictions[1:10]),
#'   train_set_correlations = cor(GROAN.KI$yield[-(1:10)], results.SNP.BLUP$predictions[-(1:10)])
#' )
#' correlations = rbind(correlations, data.frame(
#'   model = 'G-BLUP',
#'   test_set_correlations = cor(GROAN.KI$yield[1:10], results.G.BLUP$predictions[1:10]),
#'   train_set_correlations = cor(GROAN.KI$yield[-(1:10)], results.G.BLUP$predictions[-(1:10)])
#' ))
#' print(correlations)
#' }
phenoRegressor.rrBLUP = function (phenotypes, genotypes=NULL, covariances=NULL, extraCovariates = NULL, ...){
  #is rrBLUP installed?
  if (!requireNamespace("rrBLUP", quietly = TRUE)) {
    stop("rrBLUP package needed for this regressor to work. Please install it.",
         call. = FALSE)
  }

  #we need either genotypes or covariances (but not both!)
  if (!xor(is.null(genotypes), is.null(covariances))){
    stop("rrBLUP regressor requires either genotypes (SNPs) or covariances (usually kinship matrices), and not both!",
         call. = FALSE)
  }

  #which function should I call?
  if (is.null(genotypes)){
    return(phenoRegressor.rrBLUP.G(phenotypes = phenotypes, covariances = covariances, extraCovariates = extraCovariates, ...))
  }else{
    return(phenoRegressor.rrBLUP.SNP(phenotypes, genotypes, extraCovariates, ...))
  }
}

#' SNP-BLUP using rrBLUP library
#'
#' Implementation of SNP-BLUP using rrBLUP library. Not to be exported.
#'
#' @param phenotypes phenotypes, a numeric array (n x 1), missing values are predicted
#' @param genotypes SNP genotypes, one row per phenotype (n), one column per marker (m), values in 0/1/2 for
#'                  diploids or 0/1/2/...ploidy for polyploids. Can be NULL if \code{covariances} is present.
#' @param extraCovariates optional extra covariates set, one row per phenotype (n), one column per covariate (w).
#'                 If NULL no extra covariates are considered.
#' @param ... extra parameters are passed to rrBLUP::mixed.solve
#'
#' @keywords internal
#'
#' @return The function returns a list with the following fields:
#' \itemize{
#'   \item \code{predictions} : an array of (k) predicted phenotypes
#'   \item \code{hyperparams} : named vector with the following keys: Vu, Ve, beta, LL
#'   \item \code{extradata}   : list with information on trained model, coming from \code{\link[rrBLUP]{mixed.solve}}
#' }
phenoRegressor.rrBLUP.SNP = function (phenotypes, genotypes, extraCovariates = NULL, ...){
  #a handy selector to separate test and train sets
  test = is.na(phenotypes)

  #training the model
  m =  rrBLUP::mixed.solve(y = phenotypes, Z = genotypes, X = extraCovariates, ...)

  #applying the model to test set to obtain predictions and obtaining:
  #m$u     marker effects
  #m$beta  fixed effects

  #since fixed effects are optional, we must compute their contribution aside
  if (is.null(extraCovariates)){
    fi = m$beta
  }else{
    fi = as.matrix(extraCovariates) %*% m$beta
  }

  #predicting test set
  preds = as.vector(fi) + (as.matrix(genotypes) %*% m$u)[,1]

  #creating the result list
  return(list(
    predictions = preds,
    hyperparams = c(Vu=m$Vu, Ve=m$Ve, beta=m$beta, LL=m$LL),
    extradata = m
  ))
}

#' G-BLUP using rrBLUP library
#'
#' This function implements G-BLUP using rrBLUP library. Not to be exported.
#'
#' @param phenotypes phenotypes, a numeric array (n x 1), missing values are predicted
#' @param covariances square matrix (n x n) of covariances.
#' @param extraCovariates optional extra covariates set, one row per phenotype (n), one column per covariate (w).
#'                 If NULL no extra covariates are considered.
#' @param ... extra parameters are passed to rrBLUP::mixed.solve
#'
#' @keywords internal
#'
#' @return The function returns a list with the following fields:
#' \itemize{
#'   \item \code{predictions} : an array of (k) predicted phenotypes
#'   \item \code{hyperparams} : named vector with the following keys: Vu, Ve, beta, LL
#'   \item \code{extradata}   : list with information on trained model, coming from \code{\link[rrBLUP]{mixed.solve}}
#' }
phenoRegressor.rrBLUP.G = function (phenotypes, covariances, extraCovariates = NULL, ...){
  #extracting sample names
  samples = rownames(covariances)

  #we need to build a dataframe with all the required data
  df = data.frame(
    'phenotypes' = phenotypes,
    'samples' = samples
  )

  #should we add extraCovariates as fixed effects?
  if(!is.null(extraCovariates)){
    df = cbind(df, extraCovariates)
  }

  #we usually use GAUSS = FALSE to use passed kinship, but
  #the user takes precedence
  user.args = list(...)
  GAUSS.val = user.args$GAUSS
  if (is.null(GAUSS.val)){
    GAUSS.val = FALSE
  }

  #if covariances are not already a matrix, this is a good moment to change
  if (!is.matrix(covariances)){
    covariances = as.matrix(covariances)
  }

  #training the model
  m =  rrBLUP::kin.blup(
    data = df, geno = 'samples', pheno = 'phenotypes', fixed = colnames(extraCovariates),
    GAUSS = GAUSS.val, K = covariances)

  #as predictions, we use the predicted genetic values, averaged over the fixed effects
  preds = m$pred

  #creating the result list
  return(list(
    predictions = preds,
    hyperparams = c(Vg=m$Vg, Ve=m$Ve),
    extradata = m
  ))
}

#' Support Vector Regression using package e1071
#'
#' This is a wrapper around several functions from \code{e1071} package (as such, it won't work if
#' e1071 package is not installed).
#' This function implements Support Vector Regressions, meaning that the data points are projected in
#' a transformed higher dimensional space where linear regression is possible.\cr
#' \cr
#' \code{phenoRegressor.SVR} can operate in three modes: run, train and tune.\cr
#' In \strong{run} mode you need to pass the function an already tuned/trained SVR model, typically
#' obtained either directly from e1071 functions (e.g. from \link[e1071]{svm}, \link[e1071]{best.svm} and so forth)
#' or from a previous run of \code{phenoRegressor.SVR} in a different mode. The passed model is applied
#' to the passed dataset and predictions are returned.\cr
#' In \strong{train} mode a SVR model will be trained on the passed dataset using the passed hyper
#' parameters. The trained model will then be used for predictions.\cr
#' In \strong{tune} mode you need to pass one or more sets of hyperparameters. The best combination of
#' hyperparameters will be selected through crossvalidation. The best performing SVR model will be used
#' for final predictions. This mode can be very slow.\cr
#' \cr
#' There is no distinction between regular covariates (genotypes) and extra
#' covariates (fixed effects) in Support Vector Regression. If extra covariates are passed, they are
#' put together with genotypes, side by side. Same thing happens with covariances matrix. This
#' can bring to the scientifically questionable but technically correct situation of regressing
#' on a big matrix made of SNP genotypes, covariances and other covariates, all collated side by side.
#' The function makes no distinction, and it's up to the user understand what is correct in each
#' specific experiment.\cr
#'
#' @param phenotypes phenotypes, a numeric array (n x 1), missing values are predicted
#' @param genotypes SNP genotypes, one row per phenotype (n), one column per marker (m), values in 0/1/2 for
#'                  diploids or 0/1/2/...ploidy for polyploids. Can be NULL if \code{covariances} is present.
#' @param covariances square matrix (n x n) of covariances. Can be NULL if \code{genotypes} is present.
#' @param extraCovariates extra covariates set, one row per phenotype (n), one column per covariate (w).
#'                 If NULL no extra covariates are considered.
#' @param mode this parameter decides what will happen with the passed dataset
#'             \itemize{
#'             \item \code{mode = "tune"} : hyperparameters will be tuned on a grid (you may want to specify
#'                                 its values using extra params) with a call to \code{e1071::tune.svm}. Use this option if
#'                                 you have no idea about the optimal choice of hyperparameters. This mode can be very slow.
#'             \item \code{mode = "train"} : an SVR will be trained on the train dataset using the passed hyperparameters
#'                                  (if you know them). This more invokes \code{e1071::train}
#'             \item \code{mode = "run"} : you already have a tuned and trained SVR (put it into \code{tuned.model}) and
#'                                want to use it. The fastest mode.
#'             }
#' @param tuned.model a tuned and trained SVR to be used for prediction. This object is only used if
#'                    \code{mode} is equal to "run".
#' @param scale.pheno if TRUE (default) the phenotypes will be scaled and centered (before tuning or before
#'                    applying the passed tuned model).
#' @param scale.geno if TRUE the genotypes will be scaled and centered (before tuning or before
#'                    applying the passed tuned model. It is usually not a good idea, since it leads to
#'                    worse results. Defaults to FALSE.
#' @param ... all extra parameters are passed to \code{e1071::svm} or \code{e1071::tune.svm}
#'
#' @return The function returns a list with the following fields:
#' \itemize{
#'   \item \code{predictions} : an array of (n) predicted phenotypes
#'   \item \code{hyperparams} : named vector with the following keys: gamma, cost, coef0, nu, epsilon. Some
#'                              of the values may not make sense given the selected model, and will contain
#'                              default values from e1071 library.
#'   \item \code{extradata}   : depending on \code{mode} parameter, \code{extradata} will contain one of the
#'                              following:
#'                              1) a SVM object returned by e1071::tune.svm, containing both
#'                              the best performing model and the description of the training process
#'                              2) a newly trained SVR model
#'                              3) the same object passed as \code{tuned.model}
#' }
#' @family phenoRegressors
#' @seealso \link[e1071]{svm}, \link[e1071]{tune.svm}, \link[e1071]{best.svm} from e1071 package
#' @export
#' @examples \dontrun{
#' ### WARNING ###
#' #The 'tuning' part of the example can take quite some time to run,
#' #depending on the computational power.
#'
#' #using the GROAN.KI dataset, we regress on the dataset and predict the first ten phenotypes
#' phenos = GROAN.KI$yield
#' phenos[1:10] = NA
#'
#' #--------- TUNE ---------
#' #tuning the SVR on a grid of hyperparameters
#' results.tune = phenoRegressor.SVR(
#'   phenotypes = phenos,
#'   genotypes = GROAN.KI$SNPs,
#'   covariances = NULL,
#'   extraCovariates = NULL,
#'   mode = 'tune',
#'   kernel = 'linear', cost = 10^(-3:+3) #SVR-specific parameters
#' )
#'
#' #examining the predictions
#' plot(GROAN.KI$yield, results.tune$predictions,
#'      main = 'Mode = TUNING\nTrain set (black) and test set (red) regressions',
#'      xlab = 'Original phenotypes', ylab = 'Predicted phenotypes')
#' points(GROAN.KI$yield[1:10], results.tune$predictions[1:10], pch=16, col='red')
#'
#' #printing correlations
#' test.set.correlation  = cor(GROAN.KI$yield[1:10], results.tune$predictions[1:10])
#' train.set.correlation = cor(GROAN.KI$yield[-(1:10)], results.tune$predictions[-(1:10)])
#' writeLines(paste(
#'   'test-set correlation :', test.set.correlation,
#'   '\ntrain-set correlation:', train.set.correlation
#' ))
#'
#' #--------- TRAIN ---------
#' #training the SVR, hyperparameters are given
#' results.train = phenoRegressor.SVR(
#'   phenotypes = phenos,
#'   genotypes = GROAN.KI$SNPs,
#'   covariances = NULL,
#'   extraCovariates = NULL,
#'   mode = 'train',
#'   kernel = 'linear', cost = 0.01 #SVR-specific parameters
#' )
#'
#' #examining the predictions
#' plot(GROAN.KI$yield, results.train$predictions,
#'      main = 'Mode = TRAIN\nTrain set (black) and test set (red) regressions',
#'      xlab = 'Original phenotypes', ylab = 'Predicted phenotypes')
#' points(GROAN.KI$yield[1:10], results.train$predictions[1:10], pch=16, col='red')
#'
#' #printing correlations
#' test.set.correlation  = cor(GROAN.KI$yield[1:10], results.train$predictions[1:10])
#' train.set.correlation = cor(GROAN.KI$yield[-(1:10)], results.train$predictions[-(1:10)])
#' writeLines(paste(
#'   'test-set correlation :', test.set.correlation,
#'   '\ntrain-set correlation:', train.set.correlation
#' ))
#'
#' #--------- RUN ---------
#' #we recover the trained model from previous run, predictions will be exactly the same
#' results.run = phenoRegressor.SVR(
#'   phenotypes = phenos,
#'   genotypes = GROAN.KI$SNPs,
#'   covariances = NULL,
#'   extraCovariates = NULL,
#'   mode = 'run',
#'   tuned.model = results.train$extradata
#' )
#'
#' #examining the predictions
#' plot(GROAN.KI$yield, results.run$predictions,
#'      main = 'Mode = RUN\nTrain set (black) and test set (red) regressions',
#'      xlab = 'Original phenotypes', ylab = 'Predicted phenotypes')
#' points(GROAN.KI$yield[1:10], results.run$predictions[1:10], pch=16, col='red')
#'
#' #printing correlations
#' test.set.correlation  = cor(GROAN.KI$yield[1:10], results.run$predictions[1:10])
#' train.set.correlation = cor(GROAN.KI$yield[-(1:10)], results.run$predictions[-(1:10)])
#' writeLines(paste(
#'   'test-set correlation :', test.set.correlation,
#'   '\ntrain-set correlation:', train.set.correlation
#' ))
#' }
phenoRegressor.SVR = function (phenotypes, genotypes, covariances, extraCovariates, mode=c('tune', 'train', 'run'), tuned.model=NULL, scale.pheno=TRUE, scale.geno=FALSE, ...){
  #is e1071 installed?
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop("e1071 package needed for this regressor to work. Please install it.",
         call. = FALSE)
  }

  #ensuring the mode is ok
  mode = match.arg(mode)

  #ensuring mode coherency
  if (mode == 'run' & is.null(tuned.model)){
    stop("No tuned.model specified, but mode is 'run'",
         call. = FALSE)
  }

  #are scaling signal booleans?
  if (!is.boolean(scale.pheno)) stop("Parameter scale.pheno should be boolean", call. = FALSE)
  if (!is.boolean(scale.geno))  stop("Parameter scale.geno should be boolean", call. = FALSE)

  #if we get to this point, we can proceed

  #covariances and extra covariates are put together with genotypes
  if (is.null(genotypes)){
    genotypes = covariances
  }else{
    if (!is.null(covariances)){
      genotypes = cbind(genotypes, covariances)
    }
  }
  if (!is.null(extraCovariates)){
    genotypes = cbind(genotypes, extraCovariates)
  }

  #tuning, training or recycling?
  if (mode == 'run'){
    #the user passed an already tuned model, let's use it!
    model = tuned.model
    extra = tuned.model
  }else{
    #building a dataframe usable by formula interface
    df.train = data.frame(y_to_by_predicted = phenotypes, genotypes)

    #scaling (array of boolean, first position describes phenotypes, other
    #position describe genotypes)
    scale.signal = rep(scale.geno, ncol(df.train))
    scale.signal[1] = scale.pheno

    #tuning or training?
    if (mode == 'tune'){
      #for tuning, to avoid e1071 to generate warnings, we remove lines containing NAs
      train.set = !is.na(phenotypes)

      training.result = e1071::tune.svm(
        y_to_by_predicted ~ .,
        data=df.train[train.set,],
        scale=scale.signal,
        ...)
      #storing results in standard named vars
      model = training.result$best.model
    }else{
      training.result = e1071::svm(
        y_to_by_predicted ~ .,
        data=df.train,
        scale=scale.signal,
        ...)
      #storing results in standard named vars
      model = training.result
    }
    extra = training.result
  }

  #retrieving kernel name
  knames = c('linear', 'polynomial', 'radial', 'sigmoid')

  #extracting hyperparams (some assume default
  #value and are not used)
  hyp = c(
    'mode' = mode,
    'kernel' = knames[model$kernel + 1],
    'gamma' = model$gamma,
    'cost' = model$cost,
    'coef0' = model$coef0,
    'nu' = model$nu,
    'epsilon' = model$epsilon
  )

  #returning predictions and stuff
  return(list(
    predictions = predict(object=model, newdata=genotypes),
    hyperparams = hyp,
    extradata = extra
  ))
}

#' Random Forest Regression using package randomForest
#'
#' This is a wrapper around \link[randomForest]{randomForest} and related functions.
#' As such, this function will not work if randomForest package is not installed.
#' There is no distinction between regular covariates (genotypes) and extra
#' covariates (fixed effects) in random forest. If extra covariates are passed, they are
#' put together with genotypes, side by side. Same thing happens with covariances matrix. This
#' can bring to the scientifically questionable but technically correct situation of regressing
#' on a big matrix made of SNP genotypes, covariances and other covariates, all collated side by side.
#' The function makes no distinction, and it's up to the user understand what is correct in each
#' specific experiment.\cr
#' \cr
#' \strong{WARNING}: this function can be *very* slow, especially when called on thousands of SNPs.
#'
#' @param phenotypes phenotypes, a numeric array (n x 1), missing values are predicted
#' @param genotypes SNP genotypes, one row per phenotype (n), one column per marker (m), values in 0/1/2 for
#'                  diploids or 0/1/2/...ploidy for polyploids. Can be NULL if \code{covariances} is present.
#' @param covariances square matrix (n x n) of covariances. Can be NULL if \code{genotypes} is present.
#' @param extraCovariates extra covariates set, one row per phenotype (n), one column per covariate (w).
#'                 If NULL no extra covariates are considered.
#' @param ntree number of trees to grow, defaults to a fifth of the number of samples (rounded
#'              up). As per \code{randomForest} documentation, it should not be set to too
#'              small a number, to ensure that every input row gets predicted at least a few times
#' @param ... any extra parameter is passed to \code{randomForest::randomForest()}
#'
#' @return The function returns a list with the following fields:
#' \itemize{
#'   \item \code{predictions} : an array of (k) predicted phenotypes
#'   \item \code{hyperparams} : named vector with the following keys: ntree (number of grown trees)
#'                              and mtry (number of variables randomly sampled as candidates at each split)
#'   \item \code{extradata}   : the object returned by \code{randomForest::randomForest()}, containing the
#'                              full trained forest and the used parameters
#' }
#' @family phenoRegressors
#' @seealso \link[randomForest]{randomForest}
#' @export
#' @examples
#' \dontrun{
#' #using the GROAN.KI dataset, we regress on the dataset and predict the first ten phenotypes
#' phenos = GROAN.KI$yield
#' phenos[1:10]  = NA
#'
#' #calling the regressor with random forest
#' results = phenoRegressor.RFR(
#'   phenotypes = phenos,
#'   genotypes = GROAN.KI$SNPs,
#'   covariances = NULL,
#'   extraCovariates = NULL,
#'   ntree = 20,
#'   mtry = 200 #randomForest-specific parameters
#' )
#'
#' #examining the predictions
#' plot(GROAN.KI$yield, results$predictions,
#'      main = 'Train set (black) and test set (red) regressions',
#'      xlab = 'Original phenotypes', ylab = 'Predicted phenotypes')
#' points(GROAN.KI$yield[1:10], results$predictions[1:10], pch=16, col='red')
#'
#' #printing correlations
#' test.set.correlation  = cor(GROAN.KI$yield[1:10], results$predictions[1:10])
#' train.set.correlation = cor(GROAN.KI$yield[-(1:10)], results$predictions[-(1:10)])
#' writeLines(paste(
#'   'test-set correlation :', test.set.correlation,
#'   '\ntrain-set correlation:', train.set.correlation
#' ))
#' }
phenoRegressor.RFR = function (
  phenotypes, genotypes, covariances, extraCovariates,
  ntree=ceiling(length(phenotypes)/5), ...){
  #is randomForest installed?
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("randomForest package needed for this regressor to work. Please install it.",
         call. = FALSE)
  }

  #covariances and extra covariates are put together with genotypes
  if (is.null(genotypes)){
    genotypes = covariances
  }else{
    if (!is.null(covariances)){
      genotypes = cbind(genotypes, covariances)
    }
  }
  if (!is.null(extraCovariates)){
    genotypes = cbind(genotypes, extraCovariates)
  }

  #building a dataframe usable by formula interface
  df.train = data.frame(y_to_by_predicted = phenotypes, genotypes)

  #training the random forest
  RF.trained_model = randomForest::randomForest(y_to_by_predicted ~ ., data=df.train, ntree=ntree, na.action=na.omit, ...)

  #extracting hyperparams
  hyp = c(
    'ntree' = RF.trained_model$ntree,
    'mtry' = RF.trained_model$mtry
  )

  #returning predictions and stuff
  return(list(
    predictions = predict(object=RF.trained_model, newdata=genotypes),
    hyperparams = hyp,
    #in extradata we just put the trained model
    extradata = RF.trained_model
  ))
}
