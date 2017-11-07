### created for cdi grant proposal
### exploring shape of rss functions for sine model
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









ii <- 1
topN <- 20
p_grid <- rev(1/omegas)

col <- "#00000030"
wid <- 600
hei <- 350

#### SINE MODEL

## downsampled
tm <- tms[[ii]]
n <- sum(vapply(tm,nrow,c(0)))
out <- pgls(tm,periods=p_grid,BCD_flag=FALSE)
rss <- out$rss_ls
png("rss_sin_down.png",width=wid,height=hei)
par(mar=c(5,5,1,1))
plot(p_grid,rss/n,xlab="Period",ylab="RSS/n",cex.lab=1.3,cex.axis=1.3,xaxs='i',col=col)
abline(v=periods[ii],col='red',lwd=2)
abline(v=p_grid[which.min(rss)],col='orange',lwd=2)
dev.off()

## full
tm <- tms_FULL[[ii]]
n <- sum(vapply(tm,nrow,c(0)))
out <- pgls(tm,periods=p_grid,BCD_flag=FALSE)
rss <- out$rss_ls
png("rss_sin_full.png",width=wid,height=hei)
par(mar=c(5,5,1,1))
plot(p_grid,rss/n,xlab="Period",ylab="RSS/n",cex.lab=1.3,cex.axis=1.3,xaxs='i',col=col)
abline(v=periods[ii],col='red',lwd=2)
abline(v=p_grid[which.min(rss)],col='orange',lwd=2)
dev.off()


#### RRL MODEL

## down
lc <- TMtoLC(tmsc[[ii]])
n <- nrow(lc)
rss <- FitTemplate(lc,omegas,tem,NN=NN,use.errors=TRUE,use.dust=FALSE)
png("rss_rr_down.png",width=wid,height=hei)
par(mar=c(5,5,1,1))
plot(p_grid,rev(rss)/n,xlab="Period",ylab="RSS/n",cex.lab=1.3,cex.axis=1.3,xaxs='i',col=col)
abline(v=periods[ii],col='red',lwd=2)
abline(v=p_grid[which.min(rev(rss))],col='orange',lwd=2)
dev.off()


## full
lc <- TMtoLC(tmsc_FULL[[ii]])
n <- nrow(lc)
rss <- FitTemplate(lc,omegas,tem,NN=NN,use.errors=TRUE,use.dust=FALSE)
png("rss_rr_full.png",width=wid,height=hei)
par(mar=c(5,5,1,1))
plot(p_grid,rev(rss)/n,xlab="Period",ylab="RSS/n",cex.lab=1.3,cex.axis=1.3,xaxs='i',col=col)
abline(v=periods[ii],col='red',lwd=2)
abline(v=p_grid[which.min(rev(rss))],col='orange',lwd=2)
dev.off()


#### plot the light-curves
lc <- TMtoLC(tmsc[[ii]])
lc_full <- TMtoLC(tmsc_FULL[[ii]])
