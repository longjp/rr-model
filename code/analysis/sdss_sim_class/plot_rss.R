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
ebv <- extcr / tem_sdss$dust['r']
tmsc <- mapply(DustCorrect,tms,ebv,MoreArgs=list(tem=tem_sdss),SIMPLIFY=FALSE)
tmsc_FULL <- mapply(DustCorrect,tms_FULL,ebv,MoreArgs=list(tem=tem_sdss),SIMPLIFY=FALSE)

## estimate periods for both the full lc (tms_FULL) and downsampled (tms)
## using both the RRL template and (period_est) and lomb-scarge (period_est_lomb)

ii <- 1
lc <- TMtoLC(tms[[ii]])
rss <- FitTemplate(lc,omegas,tem_sdss,NN=NN,use.errors=TRUE,use.dust=TRUE)

png("freq_rss.png",width=1000,height=500)
par(mar=c(5,5,1,1))
plot(omegas,rss,xlab=expression(omega),ylab="RSS",cex.lab=2,cex.axis=2,xaxs='i')
abline(v=omegas[which.min(rss)])
dev.off()


