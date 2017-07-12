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

##################################################################
# In this file all the functionalities required to run the tests #
##################################################################

#' Compare Genomic Regressors on a Noisy Dataset
#'
#' This function runs the experiment described in a \link[=createWorkbench]{GROAN.Workbench} object on
#' the data contained in a \link[=createNoisyDataset]{GROAN.NoisyDataSet} object. It returns a \code{GROAN.Result}
#' object, which have a \link[=summary.GROAN.Result]{summary} function for quick inspection and can be
#' fed to \link{plotResult} for visual comparisons. The experiment statistics are computed via \link{measurePredictionPerformance}.\cr
#' Each time this function is invoked it will refer to a \code{runId} - an alphanumeric string identifying
#' each specific run. The \code{runId} is usually generated internally, but it is possible to pass it if
#' the intention is to join results from different runs for analysis purposes.
#'
#' @param nds a GROAN.NoisyDataSet object, containing the data (genotypes, phenotypes and so forth) plus a \code{noiseInjector} function
#' @param wb a GROAN.Workbench object, containing the regressors to be tested together with the description of the experiment
#' @param run.id an alphanumeric string identifying this specific run. If not passed it is generated using \link{createRunId}
#'
#' @return a \code{GROAN.Result} object
#' @seealso \link{measurePredictionPerformance}
#' @export
#'
#' @examples
#' \dontrun{
#' #Complete examples are found in the vignette
#' vignette('GROAN.vignette', package='GROAN')
#'
#' #Minimal example
#' #1) creating a noisy dataset with normal noise
#' nds = createNoisyDataset(
#'   name = 'PEA, normal noise',
#'   genotypes = GROAN.pea.SNPs,
#'   phenotypes = GROAN.pea.yield,
#'   noiseInjector = noiseInjector.norm,
#'   mean = 0,
#'   sd = sd(GROAN.pea.yield) * 0.5
#' )
#'
#' #2) creating a GROAN.WorkBench using default regressor and crossvalidation preset
#' wb = createWorkbench()
#'
#' #3) running the experiment
#' res = GROAN.run(nds, wb)
#'
#' #4) examining results
#' summary(res)
#' plotResult(res)
#' }
GROAN.run = function(nds, wb, run.id = createRunId()){
  #steps for user interface
  st.tot = ceiling(wb$reps/20)
  st.current = 0

  #some parameters are copied in local variables, for simplicity
  folds = wb$folds            #crossvalidation folds
  reps = wb$reps              #number of general repetitions
  strat = wb$stratified       #stratified crossvalidation?
  r.functions = wb$regressor  #regressor functions
  r.names = wb$regressor.name #regressor names
  r.parms = wb$extraParms     #regressor params
  n = length(nds$phenos)      #number of samples
  m = ncol(nds$genos)         #number of markers
  if (is.null(m)){
    m = NA
  }
  extraCovNames = getExtraCovariatesNames(nds)    #names of covariates

  #should we save the accuracies to file?
  if (is.null(wb$outfolder)){
    #no saving is required
    accuracy.save = FALSE
  }else{
    #saving
    accuracy.save = TRUE
    #let's be sure there is room for saving files
    dir.create(wb$outfolder, showWarnings = FALSE, recursive = TRUE)
    #accuracy output file name
    accuracy.filename = file.path(wb$outfolder, 'accuracy.csv')
    #flag to mark if output file already exists (and to trigger append)
    accuracy.append = file.exists(accuracy.filename)
  }

  #message to user
  writeLines(paste('GROAN starts to run, dataset', nds$name))

  #if a stratified crossvalidation is required, but no strata are present,
  #all data points will belong to dummy strata
  if(strat & is.null(nds$strata)){
    nds$strata = rep('dummyStrata', n)
  }

  #room for results
  res = NULL

  #for each required repetition
  for (i in 1:reps){
    #user interface
    st.current = st.current + 1
    if (st.current == st.tot){
      writeLines(paste(sep='', '[', run.id, '] Repetition ', i, '/', reps))
      st.current = 0
    }

    #generate the selectors for crossvalidation folds, stratified or standard
    if (strat){
      f.selectors = stratified_crossvalidation_folds(nds$strata, folds)
    }else{
      f.selectors = sample(rep(1:folds, length.out = n))
    }

    #generate noisy dataset
    noisy = getNoisyPhenotype(nds)

    #for each regressor
    for(r.curr in 1:length(r.names)){
      #room for accuracy measurement
      acc = NULL
      #for each fold
      for(f in 1:folds){
        #the samples in test set have phenotype removed (put to NA)
        phenos.train = noisy
        phenos.train[f.selectors == f] = NA

        #building the list of arguments for regressor function
        r.args = list(
          phenotypes = phenos.train,
          genotypes = nds$genos,
          covariances = nds$covs,
          extraCovariates = nds$extraCovs)
        r.args = c(r.args, r.parms[r.curr][[1]])

        #train & predict, while measuring time
        times = system.time({preds = do.call(r.functions[[r.curr]], args = r.args)})

        #ensuring that all optional return argument are there anyway
        if (is.null(preds$hyperparams)){
          preds$hyperparams = c()
        }
        if (is.null(preds$extradata)){
          preds$extradata = list()
        }

        #measure accuracy per strata, if stratified
        if (strat){
          #for each strata present in the fold
          strata.currentfold = nds$strata[f.selectors == f]
          for (s in unique(strata.currentfold)){
            #extracting the prediction for this stratum only
            preds.curr = preds$predictions[(f.selectors == f) & (nds$strata == s)]
            #extracting the original data for this stratum only
            orig.curr = nds$phenos[(f.selectors == f) & (nds$strata == s)]

            #storing the accuracies
            acc = rbind(acc, data.frame(
              'strata' = s,
              'time_per_fold' = as.numeric(times['elapsed']),
              t(measurePredictionPerformance(orig.curr, preds.curr))))
          }
        }

        #measuring accuracy regardless of strata
        acc = rbind(acc, data.frame(
          'strata' = 'no_strata',
          'time_per_fold' = as.numeric(times['elapsed']),
          t(measurePredictionPerformance(
            nds$phenos[f.selectors == f],
            preds$predictions[f.selectors == f]))))

        #if required, save to file hyperparams
        if (wb$saveHyperParms){
          #building filename/path
          hp.filename = paste(sep='', 'hyperparameters_', r.names[r.curr], '.csv')
          hp.filename = file.path(wb$outfolder, hp.filename)

          #we should put header only if file does not exist already
          newfile = !file.exists(hp.filename)

          #adding identifier for this run, a date, and number of samples and markers for
          #this
          write.table(
            x = t(c(
              run=run.id, date=format.Date(Sys.time()), repetition=i, fold=f,
              train_samples = nrow(r.args$train.geno),
              test_samples = nrow(r.args$test.geno),
              markers = m,
              extra_covariates = extraCovNames,
              preds$hyperparams)),
            file = hp.filename,
            sep = ',',
            row.names = FALSE,
            append = !newfile,
            col.names = newfile
          )
        }

        #if required, save to file extra results
        if (wb$saveExtraData){
          #let's be sure there is room for saving files
          extradata.folder = file.path(wb$outfolder, 'extra_data')
          dir.create(extradata.folder, showWarnings = FALSE, recursive = TRUE)

          #building filename
          ed.filename = paste(sep='_', r.names[r.curr], run.id, 'rep', i, 'fold', f)
          ed.filename = file.path(extradata.folder, ed.filename)

          #all the extra data go into a list
          extradata = list(preds$extradata)

          #adding real and predicted phenotypes
          extradata$GROAN_phenotypes_real = nds$phenos[f.selectors == f]
          extradata$GROAN_phenotypes_predicted = preds$predictions[f.selectors == f]

          #saving the extra data
          #save(extradata, file = ed.filename, envir = as.environment(preds))
          save(extradata, file = ed.filename)
        }
      }
      #average accuracy over folds, by strata, adding qualificators of this run
      acc.tot = data.frame(
        dataset = nds$name,
        samples = n,
        markers = m,
        extra_covariates = extraCovNames,
        regressor = r.names[r.curr],
        repetition = i,
        folds = folds,
        plyr::ddply(acc, 'strata', function(x){
          x$strata = NULL
          return(colMeans(x, na.rm = TRUE))
        })
      )

      #put in result object
      res = rbind(res, acc.tot)

      #if required, save to file accuracy
      if (accuracy.save){
        write.table(acc.tot, accuracy.filename, append=accuracy.append, col.names = !accuracy.append, sep=',', row.names = FALSE)
        accuracy.append = TRUE
      }
    }
  }

  #adding the proper class to the result to allow extra methods (summary and so forth)
  class(res) = c('GROAN.Result', class(res))

  #and we are done
  return(res)
}

#' Generate a random run id
#'
#' This function returns a partially random alphanumeric string that can be
#' used to identify a single run.
#'
#' @return a partially random alphanumeric string
#' @export
createRunId = function(){
  #generating a unique id for this run
  return(basename(tempfile(pattern='GROAN_run_', tmpdir = '')))
}

#' Summary of GROAN.Result
#'
#' Performance metrics are averaged over repetitions, so that a data.frame is produced
#' with one row per dataset/regressor/extra_covariates/strata/samples/markers/folds combination.
#'
#' @param object an object returned from \link{GROAN.run}
#' @param ... additional arguments ignored, added for compatibility to generic \code{summary} function
#'
#' @return a data.frame with averaged statistics
#' @export
summary.GROAN.Result = function(object, ...){
  #returning a summary anyway
  res.sum = plyr::ddply(.data = object, .variables = c('dataset', 'regressor', 'extra_covariates', 'strata', 'samples', 'markers', 'folds'), .fun = function(x){
    #columns to be averaged
    colset = setdiff(colnames(x), c('dataset', 'regressor', 'extra_covariates', 'strata', 'repetition'))

    return(colMeans(x[,colset]))
  })

  return (res.sum)
}
