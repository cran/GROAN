<!-- README.md is generated from README.Rmd. Please edit that file -->
GROAN
=====

GROAN provides a workbench to compare the performances of different genomic regression models. GROAN also allows to study the effect of different kinds of noise on the regression accuracies. You have to input your data (phenotypes, genotypes/covariances, other covariates...), select your regressor(s) and optionally your noise injector. GROAN is crossvalidation-oriented (but masks all the related gritty details). Output is produced in numeric form and, if package ggplot2 is installed, as plots. A small working dataset (GROAN.pea) is included and documented.

Installation
------------

GROAN is installed as a standard R package. It leverages several other packages to implement genomic regressions. These packages are not installed by default, but set as *suggested*. GROAN will ask you to install any missing package as soon as you try to use it.

Dependencies
------------

GROAN imports: plyr, rmarkdown, rrBLUP GROAN suggests: BGLR, e1071, ggplot2, knitr, randomForest

Documentation
-------------

Please see the package vignette for a complete tutorial. What follow is a minimal working example to give the gist of the tool.

``` r
#1) creating a noisy dataset with normal noise
nds = createNoisyDataset(
  name = 'PEA, normal noise',
  genotypes = GROAN.pea.SNPs,
  phenotypes = GROAN.pea.yield,
  noiseInjector = noiseInjector.norm,
  mean = 0,
  sd = sd(GROAN.pea.yield) * 0.5
)

#2) creating a GROAN.WorkBench using default regressor and crossvalidation preset
wb = createWorkbench()

#3) running the experiment
res = GROAN.run(nds, wb)

#4) examining results
summary(res)
plotResult(res)
```
