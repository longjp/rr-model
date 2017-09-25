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
##source("../params.R")

## load results from sim
load("0-fit.RData")

## relies on fact that RRL are first in list
N <- sum(cl=="rr")
N <- nrow(period_est_new)


periods <- periods[1:N]


## make 6 pairs of scatterplots of
## 1. true periods vs. old/new on down,up
## 2. sesar distances vs. old/new on down,up
## 3. schlegel dust vs. old/new on down,up



####### period estimates
lim <- c(0.2,1)
cex.lab=1.3
par(mar=c(5,5,3,0),mfcol=c(2,2))
plot(periods,period_est_new[,1],xlim=lim,ylim=lim,
     xlab="True Period",ylab=expression(bold("Downsampled Light Curve")),main="Abs. Mag-Period Dep.",cex.lab=cex.lab)
par(mar=c(5,5,0,0))
plot(periods,period_est_new_FULL[,1],xlim=lim,ylim=lim,
     xlab="True Period",ylab=expression(bold("Well Sampled Light Curve")),cex.lab=cex.lab)
par(mar=c(5,5,3,1))
plot(periods,period_est_old[,1],xlim=lim,ylim=lim,main="Fixed Absolute Mag",
     xlab="True Period",ylab="",cex.lab=cex.lab)
par(mar=c(5,5,0,1))
plot(periods,period_est_old_FULL[,1],xlim=lim,ylim=lim,
     xlab="True Period",ylab="",cex.lab=cex.lab)



## compute dust and distance estimates

DustDistance <- function(p_ests,tms,tem){
    N <- length(p_ests)
    coeffs <- matrix(0,ncol=4,nrow=N)
    for(ii in 1:N){
        tm <- tms[[ii]]
        omega <- 1/p_ests[ii]
        lc <- TMtoLC(tm)
        coes <- ComputeCoeffs(lc,omega,tem)
        coeffs[ii,] <- coes
    }
    colnames(coeffs) <- c("mu","E[B-V]","a","phi")
    d <- (10^(coeffs[,1]/5 + 1)) / 1000
    out <- cbind(d,coeffs[,2])
    colnames(out) <- c("distance","E[B-V]")
    return(out)
}

new_down <- DustDistance(period_est_new,tms,tem)
new_full <- DustDistance(period_est_new_FULL,tms_FULL,tem)
old_down <- DustDistance(period_est_old,tms,tem_old)
old_full <- DustDistance(period_est_old_FULL,tms_FULL,tem_old)



## obtain true distances and dust from sesar table



####### distance estimates






####### dust estimates
