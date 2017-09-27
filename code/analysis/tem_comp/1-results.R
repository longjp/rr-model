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


ScatterMatrix <- function(x,y1,y2,y3,y4,
                          xlab="",r1lab="",r2lab="",
                          c1lab="",c2lab="",lim=c(0.2,1)){
    cex.lab=1.3
    par(mar=c(5,5,3,0),mfcol=c(2,2))
    plot(x,y1,xlim=lim,ylim=lim,
         xlab=xlab,ylab=r1lab,main=c1lab,cex.lab=cex.lab)
    par(mar=c(5,5,0,0))
    plot(x,y2,xlim=lim,ylim=lim,
         xlab=xlab,ylab=r2lab,cex.lab=cex.lab)
    par(mar=c(5,5,3,1))
    plot(x,y3,xlim=lim,ylim=lim,main=c2lab,
         xlab=xlab,ylab="",cex.lab=cex.lab)
    par(mar=c(5,5,0,1))
    plot(x,y4,xlim=lim,ylim=lim,
         xlab=xlab,ylab="",cex.lab=cex.lab)
}

              
    

####### period estimates
ScatterMatrix(periods,period_est_new[,1],period_est_new_FULL[,1],
              period_est_old[,1],period_est_old_FULL[,1],
              xlab="True Period",r1lab=expression(bold("Downsampled Light Curve")),
              r2lab=expression(bold("Well Sampled Light Curve")),
              c1lab="Abs. Mag-Period Dep.",
              c2lab="Fixed Absolute Mag")


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
    out <- cbind(d,coeffs[,2]*tem$dust['r'])
    colnames(out) <- c("distance","E[B-V]")
    return(out)
}

new_down <- DustDistance(period_est_new[,1],tms,tem)
new_full <- DustDistance(period_est_new_FULL[,1],tms_FULL,tem)
old_down <- DustDistance(period_est_old[,1],tms,tem_old)
old_full <- DustDistance(period_est_old_FULL[,1],tms_FULL,tem_old)



## obtain true distances and dust from sesar table
sesar <- cbind(distance[1:N],extcr[1:N])


####### distance estimates
ScatterMatrix(sesar[,1],new_down[,1],new_full[,1],
              old_down[,1],old_full[,1],
              xlab="Sesar Distance",
              r1lab=expression(bold("Downsampled Light Curve")),
              r2lab=expression(bold("Well Sampled Light Curve")),
              c1lab="Abs. Mag-Period Dep.",
              c2lab="Fixed Absolute Mag",lim=c(0,120))

####### dust estimates
lim <- range(c(sesar[,2],new_down[,2],new_full[,2],old_down[,2],old_full[,2]))
ScatterMatrix(sesar[,2],new_down[,2],new_full[,2],
              old_down[,2],old_full[,2],
              xlab="Schlegel r-band Extinction",
              r1lab=expression(bold("Downsampled Light Curve")),
              r2lab=expression(bold("Well Sampled Light Curve")),
              c1lab="Abs. Mag-Period Dep.",
              c2lab="Fixed Absolute Mag",lim=lim)
