## results go here
rm(list=ls())
source("../../common/funcs.R")
source("../../fit_template/fit_template.R")
source("../funcs.R")

load("0-fit.RData")
load("../../data/clean/des.RData")

## load templates
load("../../fit_template/template_des.RData")
load("../../fit_template/template_sdss.RData")

## = period estimate analysis
## scatterplot, fraction correct, fraction correct by # epochs

## create dust corrected tms
ebv <- extcr / tem_sdss$dust['r']
tmsc_des <- mapply(DustCorrect,tms_des,ebv,MoreArgs=list(tem=tem_des),SIMPLIFY=FALSE)



pdf("period_accuracy.pdf")
lim <- c(.4,.9)
par(mar=c(5,5,1,1))
plot(periods,period_est_des[,1],xlim=lim,ylim=lim,
     xlab="True Period",ylab="DES Estimate",cex.lab=1.3)
abline(a=0,b=1)
dev.off()



head(period_est_des)

plot(density(period_est_des[,1] - period_est_des[,2],bw="SJ"))
