rm(list=ls())
library('parallel')
load("../data/clean/sdss_sim.RData")
load("../fit_template/template.RData")
source("../common/funcs.R")


## parameters for simulation
N <- 2
NN <- 10
omegas <- GetFreqs(0.2,1)
mc.cores <- 1

## estimate periods
period_est <- mclapply(1:N,FitTemplateParallel,mc.cores=mc.cores)
period_est <- matrix(unlist(period_est),ncol=5,byrow=TRUE)

save(period_est,file="estimate_params.RData")
