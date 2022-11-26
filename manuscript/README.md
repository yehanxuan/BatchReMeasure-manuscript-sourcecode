# BatchReMeasure 
Batch effect correction with sample remeasurement in completely confounded case-control 

## Overview 
BatchReMeasure implements batch effect correction with remeasured samples for a completely confounded case-control study. This repo contains the codes for our manuscript. The main functions are in the file ‘code’. All the codes to reproduce the results in the manuscript are contained in ‘Simulation’ and ‘RealData’. The users can find the scripts corresponding to the figure number in the manuscript. Our corresponding R package for this repo is called ['BatchReMeasure'](https://github.com/yehanxuan/BatchReMeasure).


## An example
I have wrapped the data generation procedure in `./Simulation/OneReplicate-New-S1.R` into the  `SimulateData` function in `./Simulation/SimFunc.R` file, which generates the data in our simulation study. Change the parameters if necessary 
```{r}
# Loading packages
source("./Simulation/SimFunc.R")
library(snowfall)
library(RConics)
library(RcppArmadillo)
library(Rcpp)
source("./code/ReMeasure_S1.R")
source("./code/Batch2Only.R")
source("./code/IgnoreBatch.R")
source("./code/Bootstrap.R")
source("./code/LSmodel.R")
source("./code/CorFunc.R")
```
Set the total sample size as 100 with 50 control and case samples. Set the standard deviation as 2; The 'Strong' correlation, 'Moderate' true effect, and 'Negative' batch effect are 0.9, 0.5, and -0.5, respectively.  
```{r}
data.obj <- SimulateData(total.size = 100, remeasure.size = 30, 
             scale = 2, correlation = 'Strong', true.effect = "Moderate", batch.effect = 'Negative')
# 0 represents control samples and 1 represents the case samples
data.obj$case
# design matrix 
data.obj$design
# response vector
data.obj$res
# index for remeasured samples 
data.obj$remeasure.ind
# response vector for remeasured samples
data.obj$remeasure.res
```
Since the batch effects and the biological effects are indistinguishable in this example, we conduct batch effect correction with remeasurement. We can compare it with the approaches that only consider the case and remeasured control samples in the second batch or ignore the batch effect. The bootstrap method like `batch.ReMeasure.S1.res` is also provided to improve the estimation of uncertainty
```{r}
remeasure.obj = batch.ReMeasure.S1(Y = data.obj$res, X = data.obj$case, Z = data.obj$design, ind.r = data.obj$remeasure.ind, Y.r = data.obj$remeasure.res )
batch2.obj = batch.Batch2.S1(Y = data.obj$res, X = data.obj$case, Z = data.obj$design, ind.r = data.obj$remeasure.ind, Y.r = data.obj$remeasure.res )
ignore.obj = batch.Ignore.S1(Y = data.obj$res, X = data.obj$case, Z = data.obj$design)
resboot.obj = batch.ReMeasure.S1.res(Y = data.obj$res, X = data.obj$case, Z = data.obj$design, ind.r = data.obj$remeasure.ind, Y.r = data.obj$remeasure.res)

remeasure.obj$a0
remeasure.obj$a0Var
remeasure.obj$a1
remeasure.obj$p.value

batch2.obj$a0
batch2.obj$a0Var
batch2.obj$p.value

ignore.obj$a0
ignore.obj$p.value
```

## Reproducibility 
The codes that generate the plots are in the file 'figure.' The users can reproduce some of the plots based on the cached data in 'data' file. Due to the storage limit, we have only uploaded data for Fig4, 5, S8, S9, S10, S16, S17, and S18. Other data and information is available on request from the author (hanxuan@tamu.edu)

Again, suppose you are able to run things on a Unix cluster. In that case, the shell scripts are already created and located in the 'Simulation' directory but do not forget to comment the parameters in the R scripts. 



