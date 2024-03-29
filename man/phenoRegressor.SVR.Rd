% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regressors.R
\name{phenoRegressor.SVR}
\alias{phenoRegressor.SVR}
\title{Support Vector Regression using package e1071}
\usage{
phenoRegressor.SVR(
  phenotypes,
  genotypes,
  covariances,
  extraCovariates,
  mode = c("tune", "train", "run"),
  tuned.model = NULL,
  scale.pheno = TRUE,
  scale.geno = FALSE,
  ...
)
}
\arguments{
\item{phenotypes}{phenotypes, a numeric array (n x 1), missing values are predicted}

\item{genotypes}{SNP genotypes, one row per phenotype (n), one column per marker (m), values in 0/1/2 for
diploids or 0/1/2/...ploidy for polyploids. Can be NULL if \code{covariances} is present.}

\item{covariances}{square matrix (n x n) of covariances. Can be NULL if \code{genotypes} is present.}

\item{extraCovariates}{extra covariates set, one row per phenotype (n), one column per covariate (w).
If NULL no extra covariates are considered.}

\item{mode}{this parameter decides what will happen with the passed dataset
\itemize{
\item \code{mode = "tune"} : hyperparameters will be tuned on a grid (you may want to specify
                    its values using extra params) with a call to \code{e1071::tune.svm}. Use this option if
                    you have no idea about the optimal choice of hyperparameters. This mode can be very slow.
\item \code{mode = "train"} : an SVR will be trained on the train dataset using the passed hyperparameters
                     (if you know them). This more invokes \code{e1071::train}
\item \code{mode = "run"} : you already have a tuned and trained SVR (put it into \code{tuned.model}) and
                   want to use it. The fastest mode.
}}

\item{tuned.model}{a tuned and trained SVR to be used for prediction. This object is only used if
\code{mode} is equal to "run".}

\item{scale.pheno}{if TRUE (default) the phenotypes will be scaled and centered (before tuning or before
applying the passed tuned model).}

\item{scale.geno}{if TRUE the genotypes will be scaled and centered (before tuning or before
applying the passed tuned model. It is usually not a good idea, since it leads to
worse results. Defaults to FALSE.}

\item{...}{all extra parameters are passed to \code{e1071::svm} or \code{e1071::tune.svm}}
}
\value{
The function returns a list with the following fields:
\itemize{
  \item \code{predictions} : an array of (n) predicted phenotypes
  \item \code{hyperparams} : named vector with the following keys: gamma, cost, coef0, nu, epsilon. Some
                             of the values may not make sense given the selected model, and will contain
                             default values from e1071 library.
  \item \code{extradata}   : depending on \code{mode} parameter, \code{extradata} will contain one of the
                             following:
                             1) a SVM object returned by e1071::tune.svm, containing both
                             the best performing model and the description of the training process
                             2) a newly trained SVR model
                             3) the same object passed as \code{tuned.model}
}
}
\description{
This is a wrapper around several functions from \code{e1071} package (as such, it won't work if
e1071 package is not installed).
This function implements Support Vector Regressions, meaning that the data points are projected in
a transformed higher dimensional space where linear regression is possible.\cr
\cr
\code{phenoRegressor.SVR} can operate in three modes: run, train and tune.\cr
In \strong{run} mode you need to pass the function an already tuned/trained SVR model, typically
obtained either directly from e1071 functions (e.g. from \link[e1071]{svm}, \link[e1071]{best.svm} and so forth)
or from a previous run of \code{phenoRegressor.SVR} in a different mode. The passed model is applied
to the passed dataset and predictions are returned.\cr
In \strong{train} mode a SVR model will be trained on the passed dataset using the passed hyper
parameters. The trained model will then be used for predictions.\cr
In \strong{tune} mode you need to pass one or more sets of hyperparameters. The best combination of
hyperparameters will be selected through crossvalidation. The best performing SVR model will be used
for final predictions. This mode can be very slow.\cr
\cr
There is no distinction between regular covariates (genotypes) and extra
covariates (fixed effects) in Support Vector Regression. If extra covariates are passed, they are
put together with genotypes, side by side. Same thing happens with covariances matrix. This
can bring to the scientifically questionable but technically correct situation of regressing
on a big matrix made of SNP genotypes, covariances and other covariates, all collated side by side.
The function makes no distinction, and it's up to the user understand what is correct in each
specific experiment.\cr
}
\examples{
\dontrun{
### WARNING ###
#The 'tuning' part of the example can take quite some time to run,
#depending on the computational power.

#using the GROAN.KI dataset, we regress on the dataset and predict the first ten phenotypes
phenos = GROAN.KI$yield
phenos[1:10] = NA

#--------- TUNE ---------
#tuning the SVR on a grid of hyperparameters
results.tune = phenoRegressor.SVR(
  phenotypes = phenos,
  genotypes = GROAN.KI$SNPs,
  covariances = NULL,
  extraCovariates = NULL,
  mode = 'tune',
  kernel = 'linear', cost = 10^(-3:+3) #SVR-specific parameters
)

#examining the predictions
plot(GROAN.KI$yield, results.tune$predictions,
     main = 'Mode = TUNING\nTrain set (black) and test set (red) regressions',
     xlab = 'Original phenotypes', ylab = 'Predicted phenotypes')
points(GROAN.KI$yield[1:10], results.tune$predictions[1:10], pch=16, col='red')

#printing correlations
test.set.correlation  = cor(GROAN.KI$yield[1:10], results.tune$predictions[1:10])
train.set.correlation = cor(GROAN.KI$yield[-(1:10)], results.tune$predictions[-(1:10)])
writeLines(paste(
  'test-set correlation :', test.set.correlation,
  '\ntrain-set correlation:', train.set.correlation
))

#--------- TRAIN ---------
#training the SVR, hyperparameters are given
results.train = phenoRegressor.SVR(
  phenotypes = phenos,
  genotypes = GROAN.KI$SNPs,
  covariances = NULL,
  extraCovariates = NULL,
  mode = 'train',
  kernel = 'linear', cost = 0.01 #SVR-specific parameters
)

#examining the predictions
plot(GROAN.KI$yield, results.train$predictions,
     main = 'Mode = TRAIN\nTrain set (black) and test set (red) regressions',
     xlab = 'Original phenotypes', ylab = 'Predicted phenotypes')
points(GROAN.KI$yield[1:10], results.train$predictions[1:10], pch=16, col='red')

#printing correlations
test.set.correlation  = cor(GROAN.KI$yield[1:10], results.train$predictions[1:10])
train.set.correlation = cor(GROAN.KI$yield[-(1:10)], results.train$predictions[-(1:10)])
writeLines(paste(
  'test-set correlation :', test.set.correlation,
  '\ntrain-set correlation:', train.set.correlation
))

#--------- RUN ---------
#we recover the trained model from previous run, predictions will be exactly the same
results.run = phenoRegressor.SVR(
  phenotypes = phenos,
  genotypes = GROAN.KI$SNPs,
  covariances = NULL,
  extraCovariates = NULL,
  mode = 'run',
  tuned.model = results.train$extradata
)

#examining the predictions
plot(GROAN.KI$yield, results.run$predictions,
     main = 'Mode = RUN\nTrain set (black) and test set (red) regressions',
     xlab = 'Original phenotypes', ylab = 'Predicted phenotypes')
points(GROAN.KI$yield[1:10], results.run$predictions[1:10], pch=16, col='red')

#printing correlations
test.set.correlation  = cor(GROAN.KI$yield[1:10], results.run$predictions[1:10])
train.set.correlation = cor(GROAN.KI$yield[-(1:10)], results.run$predictions[-(1:10)])
writeLines(paste(
  'test-set correlation :', test.set.correlation,
  '\ntrain-set correlation:', train.set.correlation
))
}
}
\seealso{
\link[e1071]{svm}, \link[e1071]{tune.svm}, \link[e1071]{best.svm} from e1071 package

Other phenoRegressors: 
\code{\link{phenoRegressor.BGLR}()},
\code{\link{phenoRegressor.RFR}()},
\code{\link{phenoRegressor.dummy}()},
\code{\link{phenoRegressor.rrBLUP}()},
\code{\link{phenoregressor.BGLR.multikinships}()}
}
\concept{phenoRegressors}
