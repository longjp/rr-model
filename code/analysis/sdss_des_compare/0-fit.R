## compare performance of old templates to new templates
## there are two comparisons
## 1. sdss (well sampled) old templates versus new templates
## 2. des (poorly sampled) old templates versus new templates
## the old templates have a slightly different shape, treat
## Y as z filter, and do not use photometric or model errors

## we run sdss/des data with old/new filter, so total of 4 runs

rm(list=ls())

## load necessary libraries
library('parallel')
source("../../fit_template/fit_template.R")
source("../../common/funcs.R")
source("../funcs.R")

## data source
load("../../data/clean/des.RData")

## parameters for simulation
tms <- rep(0,length(tms_sdss)) ## for compatibility with following source function
source("../params.R")

## run 1: sdss data with old templates
load("../../fit_template/template_sdss_old.RData")
period_est_sdss_old <- mclapply(1:N,FitTemplateParallel,
                                tms=tms_sdss,omegas=omegas,tem=tem,NN=NN,use.errors=FALSE,
                                use.dust=TRUE,topN=topN,mc.cores=mc.cores)
period_est_sdss_old <- matrix(unlist(period_est_sdss_old),ncol=topN,byrow=TRUE)

## run 2: sdss data with new templates
load("../../fit_template/template_sdss.RData")
period_est_sdss <- mclapply(1:N,FitTemplateParallel,
                            tms=tms_sdss,omegas=omegas,tem=tem,NN=NN,use.errors=TRUE,
                            use.dust=TRUE,topN=topN,mc.cores=mc.cores)
period_est_sdss <- matrix(unlist(period_est_sdss),ncol=topN,byrow=TRUE)

## run 3: des data with new templates
load("../../fit_template/template_des.RData")
period_est_des <- mclapply(1:N,FitTemplateParallel,
                            tms=tms_des,omegas=omegas,tem=tem,NN=NN,use.errors=TRUE,
                            use.dust=TRUE,topN=topN,mc.cores=mc.cores)
period_est_des <- matrix(unlist(period_est_des),ncol=topN,byrow=TRUE)


## run 4: des data with old templates
load("../../fit_template/template_sdss_old.RData")
## change Y band to z
for(ii in 1:length(tms_des)){
    tms_des[[ii]]$z <- rbind(tms_des[[ii]]$z,tms_des[[ii]]$Y)
    tms_des[[ii]]$Y <- NULL
}
period_est_des_old <- mclapply(1:N,FitTemplateParallel,
                                tms=tms_des,omegas=omegas,tem=tem,NN=NN,use.errors=FALSE,
                                use.dust=TRUE,topN=topN,mc.cores=mc.cores)
period_est_des_old <- matrix(unlist(period_est_des_old),ncol=topN,byrow=TRUE)


save(period_est_sdss_old,period_est_sdss,period_est_des_old,period_est_des,file="0-fit.RData")
