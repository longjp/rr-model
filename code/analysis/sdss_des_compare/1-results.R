## results go here
rm(list=ls())
source("../../common/funcs.R")
load("0-fit.RData")
load("../../data/clean/des.RData")

periods <- periods[1:nrow(period_est_sdss)]


lim <- c(.2,1)
par(mfrow=c(2,2))
plot(periods,period_est_des_old[,1],xlim=lim,ylim=lim,main="DES Old")
abline(a=0,b=1)
plot(periods,period_est_des[,1],xlim=lim,ylim=lim,main="DES New")
abline(a=0,b=1)
plot(periods,period_est_sdss_old[,1],xlim=lim,ylim=lim,main="Sloan Old ")
abline(a=0,b=1)
plot(periods,period_est_sdss[,1],xlim=lim,ylim=lim,main="Sloan New")
abline(a=0,b=1)



fraction_wrong <- function(est,truth,thres=.01){
    mean(abs(est-truth)/truth > .01)
}

fraction_wrong(period_est_des_old[,1],periods)
1-fraction_wrong(period_est_des[,1],periods)
fraction_wrong(period_est_sdss_old[,1],periods)
fraction_wrong(period_est_sdss[,1],periods)


lcs_des <- lapply(tms_des,TMtoLC)
n_epochs <- vapply(lcs_des,nrow,c(0))
hist(n_epochs)
?cut
cut(n_epochs,breaks=c(5,10,15,20,25,30,Inf))


a <- findInterval(n_epochs,c(10,15,20,25,30))
b <- abs(period_est_des[,1]-periods)/periods < .01
tapply(1*b,list(a),FUN=mean)
table(a)
tapply(1*b,list(a),FUN=sum)
unique(a)
?findInterval


sum(n_epochs <= 10)
sum(n_epochs < 10)
