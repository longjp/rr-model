rm(list=ls())

## load necessary libraries
library('parallel')
library('multiband')
load("../../fit_template/template_sdss.RData")
source("../../fit_template/fit_template.R")
source("../../common/funcs.R")
source("../funcs.R")

## data source
load("../../data/clean/sdss_sim_class.RData")

## parameters for simulation
source("../params.R")

## create dust corrected tms
ebv <- extcr / tem$dust['r']
tmsc <- mapply(DustCorrect,tms,ebv,MoreArgs=list(tem=tem),SIMPLIFY=FALSE)
tmsc_FULL <- mapply(DustCorrect,tms_FULL,ebv,MoreArgs=list(tem=tem),SIMPLIFY=FALSE)

## estimate periods for both the full lc (tms_FULL) and downsampled (tms)
## using both the RRL template and (period_est) and lomb-scarge (period_est_lomb)
period_est <- mclapply(1:N,FitTemplateParallel,
                       tms=tmsc,omegas=omegas,tem=tem,NN=NN,use.dust=FALSE,topN=topN,
                       mc.cores=mc.cores)
period_est <- matrix(unlist(period_est),ncol=topN,byrow=TRUE)

period_est_FULL <- mclapply(1:N,FitTemplateParallel,
                       tms=tmsc_FULL,omegas=omegas,tem=tem,NN=NN,use.dust=FALSE,topN=topN,
                       mc.cores=mc.cores)
period_est_FULL <- matrix(unlist(period_est_FULL),ncol=topN,byrow=TRUE)


period_est_lomb <- mclapply(1:N,FitLombParallel,
                            tms=tms,omegas=omegas,topN=topN,
                            mc.cores=mc.cores)
period_est_lomb <- matrix(unlist(period_est_lomb),ncol=topN,byrow=TRUE)

period_est_lomb_FULL <- mclapply(1:N,FitLombParallel,
                            tms=tms_FULL,omegas=omegas,topN=topN,
                            mc.cores=mc.cores)
period_est_lomb_FULL <- matrix(unlist(period_est_lomb_FULL),ncol=topN,byrow=TRUE)


save(period_est,period_est_lomb,period_est_FULL,period_est_lomb_FULL,file="0-fit.RData")
