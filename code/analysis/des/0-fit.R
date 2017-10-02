## performance of des templates on des RRL cross-matched with SDSS 
rm(list=ls())

## load necessary libraries
library('parallel')
source("../../fit_template/fit_template.R")
source("../../common/funcs.R")
source("../funcs.R")

## data source
load("../../data/clean/des.RData")
## templates
load("../../fit_template/template_des.RData")

## parameters for simulation
tms <- rep(0,length(tms_des)) ## for compatibility with following source function
source("../params.R")

## estimate periods with des data
omegas <- GetFreqs(0.4,0.95)
period_est_des <- mclapply(1:N,FitTemplateParallel,
                            tms=tms_des,omegas=omegas,tem=tem,NN=NN,use.errors=TRUE,
                            use.dust=TRUE,topN=topN,mc.cores=mc.cores)
period_est_des <- matrix(unlist(period_est_des),ncol=topN,byrow=TRUE)

## save output
save(period_est_des,file="0-fit.RData")
