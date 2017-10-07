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

####################### MAIN FUNCTIONS  #######################
#' Compare Genomic Regressors on a Noisy Dataset
#'
#' This function runs the experiment described in a \link[=createWorkbench]{GROAN.Workbench} object,
#' training regressor(s) on the data contained in a \link[=createNoisyDataset]{GROAN.NoisyDataSet} object
#' via parameter \code{nds}. The prediction accuracy is estimated either through crossvalidation
#' or on separate test dataset supplied via parameter \code{nds.test}.
#' It returns a \code{GROAN.Result} object, which have a \link[=summary.GROAN.Result]{summary}
#' function for quick inspection and can be fed to \link{plotResult} for visual comparisons.
#' In case of crossvalidation the test dataset in the result object will report the \code{[CV]}
#' suffix.\cr
#' The experiment statistics are computed via \link{measurePredictionPerformance}.\cr
#' Each time this function is invoked it will refer to a \code{runId} - an alphanumeric string identifying
#' each specific run. The \code{runId} is usually generated internally, but it is possible to pass it if
#' the intention is to join results from different runs for analysis purposes.
#'
#' @param nds a GROAN.NoisyDataSet object, containing the data (genotypes, phenotypes and so forth) plus a \code{noiseInjector} function
#' @param wb a GROAN.Workbench object, containing the regressors to be tested together with the description of the experiment
#' @param nds.test either a GROAN.NoisyDataSet or a list of GROAN.NoisyDataSet. The regression model(s) trained
#' on \code{nds} will be tested on \code{nds.test}
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
#'   name = 'PEA KI, normal noise',
#'   genotypes = GROAN.KI$SNPs,
#'   phenotypes = GROAN.KI$yield,
#'   noiseInjector = noiseInjector.norm,
#'   mean = 0,
#'   sd = sd(GROAN.KI$yield) * 0.5
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
GROAN.run = function(nds, wb, nds.test = NULL, run.id = createRunId()){
  #checking inputs
  check.nds.train(nds)
  check.nds.test(nds.test = nds.test, nds.train = nds)
  check.wb(wb, nds.test)

  #steps for user interface
  st.tot = ceiling(wb$reps/20)
  st.current = 0

  #some parameters are copied in local variables, for simplicity
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

  #if we don't do crossvalidation we put actual number of folds to 1
  folds.used = wb$folds
  folds.label = wb$folds
  if(is.null(folds.used)){
    folds.used = 1
    folds.label = 0
  }

  #if present, nds.test MUST be a list of GROAN.NoisyDataset objects
  if(!is.null(nds.test)){
    if ('GROAN.NoisyDataset' %in% class(nds.test)){
      nds.test = list(nds.test)
    }
  }

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

  #initial message to user
  msg.all = 'GROAN starts to run'
  msg = paste(' - training on dataset:', nds$name)
  if (folds.used == 1){
    msg = paste(sep='', msg, ' (using full dataset, NO crossvalidation)')
  }else{
    msg = paste(sep='', msg, ' (', folds.used, '-folds')
    if(strat){
      msg = paste(msg, 'stratified crossvalidation)')
    }else{
      msg = paste(msg, 'crossvalidation)')
    }
  }
  msg.all = c(msg.all, msg)

  msg = paste(' - testing on dataset(s):')
  datasets = c()
  if (folds.used != 1){
    datasets = c(datasets, paste(nds$name, '(crossvalidation)'))
  }
  if(!is.null(nds.test)){
    for(i in 1:length(nds.test)){
      datasets = c(datasets, nds.test[[i]]$name)
    }
  }
  datasets = paste(collapse = ', ', datasets)
  msg = paste(msg, datasets)
  msg.all = c(msg.all, msg)

  writeLines(msg.all)

  #room for results
  res=NULL

  #building structures for optional extra test datasets
  extra.test = list()
  if(!is.null(nds.test)){
    for(i in 1:length(nds.test)){
      extra.test$datasets  = c(extra.test$datasets,      rep(nds.test[[i]]$name, length(nds.test[[i]]$phenos)))
      extra.test$phenos    = c(extra.test$phenos,        nds.test[[i]]$phenos)
      extra.test$NA.phenos = c(extra.test$NA.phenos,     rep(NA, length(nds.test[[i]]$phenos)))
      extra.test$genos     = rbind(extra.test$genos,     nds.test[[i]]$genos)
      extra.test$covs      = rbind(extra.test$covs,      nds.test[[i]]$covs)
      extra.test$extraCovs = rbind(extra.test$extraCovs, nds.test[[i]]$extraCovs)
      extra.test$strata    = c(extra.test$strata,        nds.test[[i]]$strata)
    }
  }

  #for each required repetition
  for (i in 1:reps){
    #user interface
    st.current = st.current + 1
    if (st.current == st.tot){
      writeLines(paste(sep='', '[', run.id, '] Repetition ', i, '/', reps))
      st.current = 0
    }

    #generate the selectors for crossvalidation folds, stratified or standard,
    #slicing training set
    if (strat){
      f.selectors = stratified_crossvalidation_folds(nds$strata, folds.used)
    }else{
      f.selectors = sample(rep(1:folds.used, length.out = n))
    }

    #generate noisy dataset, for training
    noisy = getNoisyPhenotype(nds)

    #for each regressor
    for(r.curr in 1:length(r.names)){
      #room for accuracy measurement
      acc=NULL
      #for each fold
      for(f in 1:folds.used){
        #the samples in test set have phenotype removed (put to NA)
        phenos.train = noisy
        if(folds.used > 1){
          #we remove samples from train set only in case of real crossvalidation (number
          #of folds > 1). Otherwise all samples are used for training
          phenos.train[f.selectors == f] = NA
        }

        #joining current train data with extra test data
        all.data = extra.test
        all.data$datasets  = c(all.data$datasets, rep(nds$name, length(nds$phenos)))
        all.data$phenos    = c(all.data$phenos,        nds$phenos)
        all.data$NA.phenos = c(all.data$NA.phenos,     phenos.train) #this is used for training
        all.data$genos     = rbind(all.data$genos,     nds$genos)
        all.data$covs      = rbind(all.data$covs,      nds$covs)
        all.data$extraCovs = rbind(all.data$extraCovs, nds$extraCovs)
        all.data$strata    = c(all.data$strata,        nds$strata)

        #building the list of arguments for regressor function
        r.args = list(
          phenotypes      = all.data$NA.phenos,
          genotypes       = all.data$genos,
          covariances     = all.data$covs,
          extraCovariates = all.data$extraCovs)
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

        #measuring accuracy, plus other, interesting fields
        acc.current = multidataset_accuracy(all.data, preds, nds$name, do_stratified = strat)
        acc.current$time_per_fold = as.numeric(times['elapsed'])
        acc.current$nsamples.train = sum(!is.na(phenos.train))
        acc.current$regressor = r.names[r.curr]
        acc.current$markers = m
        acc.current$extra_covariates = extraCovNames
        acc.current$repetition = i
        acc.current$folds = folds.label

        acc = rbind(acc, acc.current)

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
      acc.avg = plyr::ddply(acc, c("dataset.train", "dataset.test", "strata", "regressor"), function(x){
        #numeric columns are averaged, the other are treated as character and collapsed together
        res = data.frame()
        for (n in colnames(x)){
          if(is.numeric(x[,n])){
            res[1, n] = mean(x[,n], na.rm = TRUE)
          }else{
            res[1, n] = paste(collapse=', ', unique(x[,n]))
          }
        }
        return(res)
      })

      #put in result object
      res = rbind(res, acc.avg)

      #if required, save to file accuracy
      if (accuracy.save){
        write.table(acc.avg, accuracy.filename, append=accuracy.append, col.names = !accuracy.append, sep=',', row.names = FALSE)
        accuracy.append = TRUE
      }
    }
  }

  #adding the proper class to the result to allow extra methods (summary and so forth)
  class(res) = c('GROAN.Result', class(res))

  #sorting columns for easier reading
  first.columns = c("regressor", "dataset.train", "dataset.test", "nsamples.test", "nsamples.train", "strata",
    "markers", "extra_covariates", "repetition", "folds")
  all.names = setdiff(colnames(res), first.columns)
  all.names = c(first.columns, all.names)
  res = res[,all.names]

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
  res.sum = plyr::ddply(.data = object, .variables = c('dataset.train', 'dataset.test', 'regressor', 'extra_covariates', 'strata', 'nsamples.train', 'markers', 'folds'), .fun = function(x){
    #columns to be averaged (removing everything else)
    colset = setdiff(colnames(x), c(
      'dataset.train', 'dataset.test',
      'regressor', 'extra_covariates', 'strata', 'repetition'))

    return(colMeans(x[,colset]))
  })

  return (res.sum)
}

# INPUT CHECKERS ----------------------------------------------------------
#fails if nds is not a GROAN.NoisyDataset object
check.nds.train = function(nds){
  if (!'GROAN.NoisyDataset' %in% class(nds)){
    stop('Parameter "nds" is not of class GROAN.NoisyDataset', call. = FALSE)
  }
}

#nds test must be a single GROAN.NoisyDataset or a list of GROAN.NoisyDataset. Moreover,
#all elements of nds.test must pass the are.compatible() test with nds.train. Finally, no
#duplicate dataset names are allowed
check.nds.test = function(nds.test, nds.train){
  #NULL is OK
  if(is.null(nds.test)){
    return()
  }

  #if it's a single GROAN.NoisyDataset we transform it into a single slotted list
  if ('GROAN.NoisyDataset' %in% class(nds.test)){
    nds.test = list(nds.test)
  }

  #now we should have, in any case, a list of GROAN.NoisyDataset.
  if (!'list' %in% class(nds.test)){
    stop('Parameter "nds.test" should be a single GROAN.NoisyDataset or a list of them', call. = FALSE)
  }

  ds.names = nds.train$name
  for(i in 1:length(nds.test)){
    #check if it's a GROAN.NoisyDataset
    if (!'GROAN.NoisyDataset' %in% class(nds.test[[i]])){
      stop('Parameter "nds.test" should be a single GROAN.NoisyDataset or a list of them', call. = FALSE)
    }

    #check if it's compatible
    if (!are.compatible(nds.test[[i]], nds.train)){
      msg = paste('Dataset', nds.test[[i]]$name, 'and', nds.train$name, 'are not compatible. See function are.compatible for details.')
      stop(msg, call. = FALSE)
    }

    #check if the name is a duplicate
    if(nds.test[[i]]$name %in% ds.names){
      stop(paste('Duplicate name of dataset:', nds.test[[i]]$name), call. = FALSE)
    }else{
      ds.names = c(ds.names, nds.test[[i]]$name)
    }
  }
}

#fails if wb is not a GROAN.Workbench object, and if no crossvalidation without test set
check.wb = function(wb, nds.test){
  if (!'GROAN.Workbench' %in% class(wb)){
    stop('Parameter "wb" is not of class GROAN.Workbench', call. = FALSE)
  }
  if (is.null(wb$folds) & is.null(nds.test)){
    stop('No test set is possible: crossvalidation is not selected and nds.test is NULL', call. = FALSE)
  }
}

# SUPPORT FUNCTIONS -------------------------------------------------------
multidataset_accuracy = function(all.data, predictions, train.dataset.name, do_stratified = FALSE){
  res = NULL

  #for each available dataset
  for (dataset.current in unique(all.data$datasets)){
    #this selector tells what samples have been predicted
    prediction.selector = all.data$datasets == dataset.current

    #if current dataset is the train one special attention to crossvalidation is required
    if(dataset.current == train.dataset.name){
      #do we have any missing data?
      prediction.selector = (all.data$datasets == dataset.current) & is.na(all.data$NA.phenos)
      #in the eventuality of returning accuracy for the current
      #dataset, we keep notes on the fact that it is the product
      #of crossvalidation
      dataset.current = paste(dataset.current, '[CV]')
    }

    #do we have something to work on?
    if(!any(prediction.selector)){
      #nothing to do here, we simply skip
      next
    }

    #at this point we can measure the accuracies
    acc = measurePredictionPerformance(
      all.data$phenos[prediction.selector],
      predictions$predictions[prediction.selector])

    res = rbind(res, data.frame(
      dataset.train = train.dataset.name,
      dataset.test = dataset.current,
      nsamples.test = sum(prediction.selector),
      strata = 'no_strata',
      t(acc)
    ))

    #doing the same thing, by strata
    if(do_stratified){
      #for each strata in the current dataset
      for(str in unique(all.data$strata[prediction.selector])){
        prediction.selector.str = prediction.selector & (all.data$strata == str)

        acc = measurePredictionPerformance(
          all.data$phenos[prediction.selector.str],
          predictions$predictions[prediction.selector.str])

        res = rbind(res, data.frame(
          dataset.train = train.dataset.name,
          dataset.test = dataset.current,
          nsamples.test = sum(prediction.selector.str),
          strata = str,
          t(acc)
        ))

      }
    }
  }

  return(res)
}
