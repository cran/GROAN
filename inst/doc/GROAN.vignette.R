## ---- eval=TRUE, include=TRUE, echo=FALSE--------------------------------
#just to be sure
library(GROAN)

## ---- eval=FALSE, include=TRUE, echo=TRUE, results='hide'----------------
#  library(GROAN)
#  
#  #array of phenotypes
#  GROAN.pea.yield
#  
#  #dataframe of SNP genotypes
#  GROAN.pea.SNPs
#  
#  #dataframe of realized genotypic kinship
#  GROAN.pea.kinship

## ---- eval=TRUE, include=TRUE, echo=TRUE---------------------------------

#creating a NoisyDataset without any extra noise injected
nds.no_noise = createNoisyDataset(
  name = 'PEA, no noise',
  genotypes = GROAN.pea.SNPs, 
  phenotypes = GROAN.pea.yield
)

#creating a NoisyDataset adding noise sampled from a normal distribution
nds.normal_noise = createNoisyDataset(
  name = 'PEA, normal noise',
  genotypes = GROAN.pea.SNPs, 
  phenotypes = GROAN.pea.yield,
  noiseInjector = noiseInjector.norm,
  mean = 0,
  sd = sd(GROAN.pea.yield) * 0.5
)

## ----eval=TRUE, echo=TRUE, fig.align="center", fig.height=5, fig.width=5, include=TRUE----
#plotting the original phenotypes
plot(GROAN.pea.yield, pch=20, main = 'True (black) vs. Noisy (red)', xlab = 'Samples', ylab = 'Phenotypes')
#plotting an instance of the phenotypes with noise injected 
points(getNoisyPhenotype(nds.normal_noise), col='red')

## ----eval=TRUE, include=TRUE, echo=TRUE, , results='hold'----------------
#average correlation oscillates around 0.89
cor(GROAN.pea.yield, getNoisyPhenotype(nds.normal_noise))
cor(GROAN.pea.yield, getNoisyPhenotype(nds.normal_noise))
cor(GROAN.pea.yield, getNoisyPhenotype(nds.normal_noise))

## ----eval=TRUE, include=TRUE, echo=TRUE----------------------------------
#no noise injector ==> the original phenotypes are returned
all(GROAN.pea.yield == getNoisyPhenotype(nds.no_noise))

## ---- eval=TRUE, include=TRUE, echo=TRUE---------------------------------
#creating a Workbench with default values 
wb = createWorkbench()

## ---- eval=TRUE, include=TRUE, echo=TRUE---------------------------------
#creating a Workbench with default values explicitly assigned 
wb = createWorkbench(
  #parameters defining crossvalidation
  folds = 5, reps = 10, stratified = FALSE, 
  
  #parameters defining save-on-hard-disk policy
  outfolder = NULL, saveHyperParms = FALSE, saveExtraData = FALSE,
  
  #a regressor
  regressor = phenoRegressor.rrBLUP, regressor.name = 'rrBLUP'
)

## ---- eval=TRUE, include=TRUE, echo=TRUE---------------------------------
#adding a regressor to an existing Workbench
wb = addRegressor(
  #the Workbench to be updater
  wb,
  #the new regressor
  regressor = phenoRegressor.BGLR, regressor.name = 'Bayesian Lasso',
  
  #regressor-specific parameters
  type = 'BL'
)

## ---- eval=FALSE, include=TRUE, echo=TRUE--------------------------------
#  #executing two GROAN test, same workbench, different data set
#  res.no_noise     = GROAN.run(nds.no_noise, wb)
#  res.normal_noise = GROAN.run(nds.normal_noise, wb)
#  
#  #putting the results together for further analysis
#  res.total = rbind(res.no_noise, res.normal_noise)

## ---- eval=FALSE, include=TRUE, echo=TRUE--------------------------------
#  #defaults is a boxplot of Pearson's correlations
#  p = plotResult(res.total)
#  print(p)

## ---- out.width = "400px", eval=TRUE, include=TRUE, echo=FALSE-----------
knitr::include_graphics('plot1.png')

## ---- eval=FALSE, include=TRUE, echo=TRUE--------------------------------
#  #a barplot with 95% confidence interval of Pearson's correlations
#  p = plotResult(res.total, plot.type = 'bar_conf95')
#  print(p)

## ---- out.width = "400px", eval=TRUE, include=TRUE, echo=FALSE-----------
knitr::include_graphics('plot2.png')

## ---- eval=FALSE, include=TRUE, echo=TRUE--------------------------------
#  #a barplot of execution times per fold, in seconds
#  p = plotResult(res.total, plot.type = 'bar', variable = 'time')
#  print(p)

## ---- out.width = "400px", eval=TRUE, include=TRUE, echo=FALSE-----------
knitr::include_graphics('plot3.png')

## ---- eval=FALSE, include=TRUE, echo=TRUE--------------------------------
#  #A noise injector adding a fixed bias to a random subset of about 50% of the data
#  my.noiseInjector.bias = function(phenotypes){
#    #an array of random booleans
#    labels = runif(length(phenotypes)) > 0.5
#  
#    #adding the bias (a fixed value of 5)
#    phenotypes[labels] = phenotypes[labels] + 5
#  
#    #returning the original phenotypes with added noise
#    return(phenotypes)
#  }

## ---- eval=FALSE, include=TRUE, echo=TRUE--------------------------------
#  #A NoisyDataSet that embeds the bias noise
#  nds.bias_noise = createNoisyDataset(
#    name = 'PEA, bias noise',
#    genotypes = GROAN.pea.SNPs,
#    phenotypes = GROAN.pea.yield,
#    noiseInjector = my.noiseInjector.bias   #the function we defined above
#  )

## ---- eval=FALSE, include=TRUE, echo=TRUE--------------------------------
#  #An improved version of the above function, this time the bias is not fixed
#  my.noiseInjector.bias2 = function(phenotypes, bias = 5){
#    #an array of random booleans
#    labels = runif(length(phenotypes)) > 0.5
#  
#    #adding the bias (from the function argument)
#    phenotypes[labels] = phenotypes[labels] + bias
#  
#    #returning the original phenotypes with added noise
#    return(phenotypes)
#  }

## ---- eval=FALSE, include=TRUE, echo=TRUE--------------------------------
#  #A NoisyDataSet with bias noise function, using the second version of the function
#  nds.bias_noise2 = createNoisyDataset(
#    name = 'PEA, bias noise, second function',
#    genotypes = GROAN.pea.SNPs,
#    phenotypes = GROAN.pea.yield,
#    noiseInjector = my.noiseInjector.bias2,   #the new version
#    bias = 20 #if omitted the default would be used
#  )

## ---- eval=FALSE, include=TRUE, echo=TRUE--------------------------------
#  phenoRegressor.dummy = function (phenotypes, genotypes, covariances, extraCovariates){
#    #no operation is actually required
#    return(list(
#      predictions = runif(length(phenotypes)), #predictions
#      hyperparams = c(some = 1, params='a'),
#      extradata = 'filler extra data for phenoRegressor.dummy'
#    ))
#  }

