## results go here
rm(list=ls())
source("../../common/funcs.R")
source("../../fit_template/fit_template.R")

load("0-fit.RData")
load("../../data/clean/des.RData")
load("../../fit_template/template_des.RData")

## = period estimate analysis
## scatterplot, fraction correct, fraction correct by # epochs

lim <- c(.4,.95)
par(mar=c(5,5,1,1))
plot(periods,period_est_des[,1],xlim=lim,ylim=lim,
     xlab="True Period",ylab="DES Estimate",cex.lab=1.3)
abline(a=0,b=1)
     
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
?findInterval


sum(n_epochs <= 15)
sum(n_epochs < 15)



## = compute coefficients using des data 
## = compare to sesar distance, schlegel ebv
## make all plots together

coeffs_des <- matrix(0,nrow=length(tms_des),ncol=4)
for(ii in 1:length(tms_des)){
    omega <- 1 / period_est_des[ii,1]
    lc <- TMtoLC(tms_des[[ii]])
    coeffs_des[ii,] <- ComputeCoeffs(lc,omega,tem)
}
colnames(coeffs_des) <- c("mu","ebv","a","phi")


d_des <- 10^(coeffs_des[,1]/5 + 1)/1000
par(mar=c(5,5,1,1))
plot(distance,d_des,xlab="Distance Sesar",ylab="Distance DES",
     cex.lab=1.3)
abline(a=0,b=1)
abline(a=0,b=1.05)
abline(a=0,b=0.95)



e_des <- coeffs_des[,2]*tem$dust['r']
par(mar=c(5,5,1,1))
lim <- c(0,max(c(extcr,e_des)))
plot(extcr,e_des,xlab="Extinction r Schlegel",ylab="DES Extinction r",
     cex.lab=1.3,xlim=lim,ylim=lim)
abline(a=0,b=1)









## compare parameter estimates to sdss well sampled
## scatter plots of mu versus mu, phase versus phase, amp versus amp
load("../../fit_template/template_sdss.RData") ## load sdss templates
coeffs_sdss <- matrix(0,nrow=length(tms_sdss),ncol=4)
for(ii in 1:length(tms_sdss)){
    omega <- 1 / periods[ii]
    lc <- TMtoLC(tms_sdss[[ii]])
    coeffs_sdss[ii,] <- ComputeCoeffs(lc,omega,tem)
}
colnames(coeffs_sdss) <- c("mu","ebv","a","phi")


ii <- 0

ii <- ii + 1
lim <- range(c(coeffs_sdss[,ii],coeffs_des[,ii]))
plot(coeffs_sdss[,ii],coeffs_des[,ii],xlim=lim,ylim=lim)
abline(a=0,b=1)


### things look pretty good except for phase, which is totally screwed up

