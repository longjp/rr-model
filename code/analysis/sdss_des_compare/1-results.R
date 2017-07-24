## results go here
rm(list=ls())
load("results.RData")
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

