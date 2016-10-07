rm(list=ls())

## load necessary libraries
library('parallel')
library('multiband')
load("../../fit_template/template.RData")
source("../../fit_template/template.R")
source("../../common/funcs.R")
source("../funcs.R")

## data source
load("../../data/clean/sdss_class.RData")

## parameters for simulation
source("../params.R")

## estimate periods
period_est <- mclapply(1:N,FitTemplateParallel,
                       tms=tms,omegas=omegas,tem=tem,NN=NN,topN=topN,
                       mc.cores=mc.cores)
period_est <- matrix(unlist(period_est),ncol=topN,byrow=TRUE)

save(period_est,file="results.RData")
