rm(list=ls())

## load necessary libraries
library('parallel')
library('multiband')
load("../fit_template/template.RData")
source("../fit_template/template.R")
source("../common/funcs.R")
source("funcs.R")

## data source
load("../data/clean/sdss_sim.RData")


## parameters for simulation
N <- 2
NN <- 10
omegas <- GetFreqs(0.2,1)
mc.cores <- 1


## TODO: parameter for number of periods to return


## estimate periods
period_est <- mclapply(1:N,FitTemplateParallel,
                       tms=tms,omegas=omegas,tem=tem,NN=NN,
                       mc.cores=mc.cores)
period_est <- matrix(unlist(period_est),ncol=5,byrow=TRUE)


period_est_lomb <- mclapply(1:N,FitLombParallel,
                            tms=tms,omegas=omegas,
                            mc.cores=mc.cores)
period_est_lomb <- matrix(unlist(period_est_lomb),ncol=5,byrow=TRUE)

save(period_est,period_est_lomb,file="sdss_sim_period_est.RData")
