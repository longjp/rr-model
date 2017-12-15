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

lim <- c(.2,1)
pdf("period_accuracy.pdf")
par(mar=c(5,5,1,1))
plot(periods,period_est_des[,1],xlim=lim,ylim=lim,
     xlab="True Period",ylab="DES Estimate",cex.lab=1.3)
abline(a=0,b=1)
dev.off()
     
fraction_wrong <- function(est,truth,thres=.01){
    mean(abs(est-truth)/truth > thres)
}

1-fraction_wrong(period_est_des[,1],periods)


lcs_des <- lapply(tms_des,TMtoLC)
n_epochs <- vapply(lcs_des,nrow,c(0))
hist(n_epochs)
?cut
cut(n_epochs,breaks=c(5,10,15,20,25,30,Inf))


a <- findInterval(n_epochs,c(9,15,20,25,30))
b <- abs(period_est_des[,1]-periods)/periods < .01
tapply(1*b,list(a),FUN=mean)
table(a)
unique(a)


sum(n_epochs <= 15)
sum(n_epochs < 15)



## = compute coefficients using des data 
## = compare to sesar distance, schlegel ebv
## make all plots together

coeffs_des <- matrix(0,nrow=length(tms_des),ncol=4)
for(ii in 1:length(tms_des)){
    omega <- 1 / period_est_des[ii,1]
    lc <- TMtoLC(tmsc_des[[ii]])
    coeffs_des[ii,] <- ComputeCoeffs(lc,omega,tem_des,use.dust=FALSE)
}
colnames(coeffs_des) <- c("mu","ebv","a","phi")


cols <- 1*(abs(period_est_des[,1] - periods) / periods > .01) + 1

####### compare distances
d_des <- 10^(coeffs_des[,1]/5 + 1)/1000
lim <- range(c(distance,d_des))
pdf("distance_des.pdf")
par(mar=c(5,5,1,1))
plot(distance,d_des,xlab="Distance Sesar",ylab="Distance DES",
     cex.lab=1.3,ylim=lim,xlim=lim,col=cols,pch=cols)
abline(a=0,b=1)
abline(a=0,b=1.05,lty=2)
abline(a=0,b=0.95,lty=2)
legend("bottomright",c("Identity","5% Scatter"),lty=1:2,lwd=2)
dev.off()

pdf("distance_des_hist.pdf")
hist(d_des/distance,xlab="DES Distance / Sesar Distance")
dev.off()
