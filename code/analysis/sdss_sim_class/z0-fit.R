## contains an alternative pipeline that
## studies feature error distributions
rm(list=ls())

## load necessary libraries
library('parallel')
load("../../fit_template/template_sdss.RData")
source("../../fit_template/fit_template.R")
source("../../common/funcs.R")
source("../funcs.R")

## data source
load("../../data/clean/sdss_sim_class.RData")

## parameters for simulation
source("../params.R")

## estimate periods for both the full lc (tms_FULL) and downsampled (tms)
## using both the RRL template
period_est <- mclapply(1:N,FitTemplateParallel,
                       tms=tms,omegas=omegas,tem=tem_sdss,NN=NN,use.dust=TRUE,topN=topN,
                       mc.cores=mc.cores)
period_est <- matrix(unlist(period_est),ncol=topN,byrow=TRUE)

period_est_FULL <- mclapply(1:N,FitTemplateParallel,
                       tms=tms_FULL,omegas=omegas,tem=tem_sdss,NN=NN,use.dust=TRUE,topN=topN,
                       mc.cores=mc.cores)
period_est_FULL <- matrix(unlist(period_est_FULL),ncol=topN,byrow=TRUE)

save(period_est,period_est_FULL,file="z0-fit.RData")
