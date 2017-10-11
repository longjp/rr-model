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

## relies on fact that RRL are first in list
N <- sum(cl=="rr")
omegas <- GetFreqs(0.4,0.95) ## frequency grid

## create dust corrected tms
ebv <- extcr / tem$dust['r']
tmsc <- mapply(DustCorrect,tms,ebv,MoreArgs=list(tem=tem),SIMPLIFY=FALSE)
tmsc_FULL <- mapply(DustCorrect,tms_FULL,ebv,MoreArgs=list(tem=tem),SIMPLIFY=FALSE)

## estimate periods for both the full lc (tms_FULL) and downsampled (tms)
## using old and new templates
period_est_new <- mclapply(1:N,FitTemplateParallel,
                           tms=tmsc,omegas=omegas,tem=tem,NN=NN,use.dust=FALSE,topN=topN,
                           mc.cores=mc.cores)
period_est_new <- matrix(unlist(period_est_new),ncol=topN,byrow=TRUE)

period_est_new_FULL <- mclapply(1:N,FitTemplateParallel,
                                tms=tmsc_FULL,omegas=omegas,tem=tem,NN=NN,use.dust=FALSE,topN=topN,
                                mc.cores=mc.cores)
period_est_new_FULL <- matrix(unlist(period_est_new_FULL),ncol=topN,byrow=TRUE)


period_est_old <-  mclapply(1:N,FitTemplateParallel,
                            tms=tms,omegas=omegas,tem=tem_old,NN=NN,topN=topN,
                            mc.cores=mc.cores)
period_est_old <- matrix(unlist(period_est_old),ncol=topN,byrow=TRUE)

period_est_old_FULL <- mclapply(1:N,FitTemplateParallel,
                                tms=tms_FULL,omegas=omegas,tem=tem_old,NN=NN,topN=topN,
                                mc.cores=mc.cores)
period_est_old_FULL <- matrix(unlist(period_est_old_FULL),ncol=topN,byrow=TRUE)


save(period_est_new,period_est_new_FULL,
     period_est_old,period_est_old_FULL,file="0-fit.RData")
