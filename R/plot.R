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

#in this file functions for graphical rendering

#' Plot results of a run
#'
#' This function uses ggplot2 package (which must be installed) to
#' graphically render the result of a run. The function receive as
#' input the output of GROAN.run and returns a ggplot2 object (that
#' can be further customized).
#' Currently implemented types of plot are:
#' \itemize{
#'   \item \code{box} : boxplot, showing the distribution of repetitions. See \link[ggplot2]{geom_boxplot}
#'   \item \code{bar} : barplot, showing the average over repetitions. See \link[ggplot2]{stat_summary}
#'   \item \code{bar_conf95} : same as 'bar', but with 95\% confidence intervals
#'   }
#'
#' @param res a result data frame containing the output of GROAN.run
#' @param variable name of the variable to be used as y values
#' @param x.label select what to put on x-axis between both train and test dataset (default), train dataset only or test dataset only
#' @param plot.type a string indicating the type of plot to be obtained
#' @param strata string determining behaviour toward strata. If \code{'no_strata'} will plot
#'               accuracies not considering strata. If \code{'avg_strata'} will average single
#'               strata accuracies. If \code{'single'} each strata will be represented separately.
#'
#' @return a ggplot2 object
#' @export
plotResult = function (res,
                       variable=c('pearson', 'spearman', 'rmse', 'time_per_fold', 'coeff_det', 'mae'),
                       x.label = c('both', 'train_only', 'test_only'),
                       plot.type=c('box', 'bar', 'bar_conf95'),
                       strata = c('no_strata', 'avg_strata', 'single')
                       ){
  #ensuring arguments consistency
  variable = match.arg(variable)
  x.label = match.arg(x.label)
  plot.type = match.arg(plot.type)
  strata = match.arg(strata)

  #is ggplot2 installed?
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package needed for this function to work. Please install it.",
         call. = FALSE)
  }

  #depending on strata we have different data manipulation
  if(strata == 'no_strata'){
    #disregarding strata
    res = subset(res, strata == 'no_strata')
  }
  if(strata == 'avg_strata'){
    #removing all strata, creating average strata res
    res = subset(res, strata != 'no_strata')
    res = plyr::ddply(
      res,
      c("dataset.train", "dataset.test", "nsamples.train", "nsamples.test", "markers", "extra_covariates", "regressor", "repetition", "folds"),
      function(x){  colMeans(x[,c("time_per_fold", "pearson", "spearman",
                             "cor_success", "rmse", "mae", "coeff_det")], na.rm = TRUE)
      })
  }
  if(strata == 'single'){
    #removing no_strata, then adding strata to dataset, so they differentiate
    res = subset(res, strata != 'no_strata')
    res$dataset.test = paste(res$dataset.test, '\nSTRATUM:', res$strata)
  }

  #creating a new column to be used for X
  if (x.label == 'both'){
    res$X = paste(sep='', 'TRAIN: ', res$dataset.train, '\nTEST: ', res$dataset.test)
  }
  if (x.label == 'train_only'){
    res$X = res$dataset.train
  }
  if (x.label == 'test_only'){
    res$X = res$dataset.test
  }

  #creating the base plot

  #adding the type
  if (plot.type == 'box'){
    #base plot object plus boxplot
    p = ggplot2::ggplot(res, ggplot2::aes(x=res$X, fill=res$regressor, y=res[,variable])) +
        ggplot2::geom_boxplot()
  }else{
    #as first step, we compute averages and confidence intervals
    lims = getConfLimits(res, variable)

    #defining the dodge
    dodge = ggplot2::position_dodge(width=0.9)

    #building the bar plot
    p = ggplot2::ggplot(lims, ggplot2::aes(x=lims$X, fill=lims$regressor, y=lims$m)) +
      ggplot2::geom_bar(position=dodge, stat="identity") +
      ggplot2::ylab(variable)

    #if necessary, we add the confidence intervals
    if (plot.type == 'bar_conf95'){
      p = p + ggplot2::geom_errorbar(ggplot2::aes(ymax=lims$ubound, ymin=lims$lbound), position=dodge, width=0.25)
    }
  }

  #adding better axis and legend labels
  p = p + ggplot2::ylab(variable) +
    ggplot2::xlab('') +
    ggplot2::guides(fill=ggplot2::guide_legend(title='Model'))

  return(p)
}

# Extract averages and confidence intervals from GROAN result
#
# Internal function. Given a res from GROAN.run() and
# the variable of interest, creates a new dataframe containing, for
# each combination of dataset and regressor, average value and
# confidence intervals at 95% under normal distribution hypothesis.
#
#' @keywords internal
#'
# @param res returned from GROAN.run()
# @param variable name of the variable to be analized
#
# @return a data frame with: 'regressor', 'X', 'm', 'lbound', 'ubound'
getConfLimits = function(res, variable){
  lims = data.frame()
  #for each combination of datasets and regressors
  for (r in unique(res$regressor)){
    for(d in unique(res$X)){
      #isolating the interesting data
      tmp = subset(res, res$regressor == r & res$X == d)[,variable]
      #mean and confidence intervals
      tmp.mean = mean(tmp)
      tmp.sd = sd(tmp)
      tmp.e = qnorm(0.975) * tmp.sd / sqrt(length(tmp))
      lims = rbind(lims, data.frame(
        regressor = r,
        X = d,
        m = tmp.mean,
        lbound = tmp.mean - tmp.e,
        ubound = tmp.mean + tmp.e
      ))
    }
  }

  #making sure we are not losing the order (in case of factors)
  if (is.factor(res$X)){
    lims$X = factor(x = lims$X, levels = levels(res$X))
  }

  return (lims)
}
