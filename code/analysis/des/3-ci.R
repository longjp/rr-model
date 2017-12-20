## create confidence sets for whether best fitting period is actually correct
## idea: 1) select p best fit frequencies
##       2) bootstrap sample light curve B times, calculate RSS at these p frequencies
##       3) calculate fraction of times omega_hat has min RSS in bootstrap
##       4) this is a confidence that the estimated period is correct
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



lim <- c(.4,.9)
par(mar=c(5,5,1,1))
plot(periods,period_est_des[,1],xlim=lim,ylim=lim,
     xlab="True Period",ylab="DES Estimate",cex.lab=1.3)
abline(a=0,b=1)





boot <- function(lc){
    return(lc[sample(1:nrow(lc),replace=TRUE),])
}


B <- 100
p <- ncol(period_est_des)
probC <- rep(0,length(tmsc_des))
for(ii in 1:length(tmsc_des)){
    print(ii)
    lc <- TMtoLC(tmsc_des[[ii]])
    omega_ests <- 1/period_est_des[ii,]
    boot_rss <- rep(0,B)
    for(jj in 1:B){
        lcb <- boot(lc)
        r <- vapply(omega_ests,function(x){min(ComputeRSSPhase(lcb,x,tem_des,use.dust=FALSE))},c(0))
        boot_rss[jj] <- r[1] - min(r[2:p])
    }
    probC[ii] <- mean(boot_rss < 0)
}



## determine if period correct, plot against correct
correct <- 1*(abs((periods-period_est_des[,1])) < 0.01)
pdf("cor_vs_prob.pdf")
plot(probC,jitter(correct),xlab="prob correct",ylab="correct (with jitter)")
dev.off()
