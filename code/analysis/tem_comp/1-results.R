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


## make 3 pairs of scatterplots of
## 1. true periods vs. old/new on down,up
## 2. sesar distances vs. old/new on down,up
## 3. schlegel dust vs. old/new on down,up
ScatterMatrix <- function(x,y1,y2,y3,y4,
                          xlab="",r1lab="",r2lab="",
                          c1lab="",c2lab="",lim=c(0.2,1),
                          y1col=1,y2col=1,y3col=1,y4col=1){
    cex.lab=1.3
    par(mar=c(5,5,3,0),mfcol=c(2,2))
    plot(x,y1,xlim=lim,ylim=lim,
         xlab=xlab,ylab=r1lab,main=c1lab,cex.lab=cex.lab,col=y1col)
    abline(a=0,b=1)
    par(mar=c(5,5,0,0))
    plot(x,y2,xlim=lim,ylim=lim,
         xlab=xlab,ylab=r2lab,cex.lab=cex.lab,col=y2col)
    abline(a=0,b=1)
    par(mar=c(5,5,3,1))
    plot(x,y3,xlim=lim,ylim=lim,main=c2lab,
         xlab=xlab,ylab="",cex.lab=cex.lab,col=y3col)
    abline(a=0,b=1)
    par(mar=c(5,5,0,1))
    plot(x,y4,xlim=lim,ylim=lim,
         xlab=xlab,ylab="",cex.lab=cex.lab,col=y4col)
    abline(a=0,b=1)
}

              
    

####### period estimates
pdf("figs/1-results-periods.pdf",height=9,width=9)
ScatterMatrix(periods,period_est_new[,1],period_est_new_FULL[,1],
              period_est_old[,1],period_est_old_FULL[,1],
              xlab="True Period",r1lab=expression(bold("Downsampled Light Curve")),
              r2lab=expression(bold("Well Sampled Light Curve")),
              c1lab="New Beta0",
              c2lab="Old Beta0")
dev.off()

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
    colnames(out) <- c("distance","extcr")
    return(out)
}


## compute dust, distance using each method / data set
new_down <- DustDistance(period_est_new[,1],tms,tem)
new_full <- DustDistance(period_est_new_FULL[,1],tms_FULL,tem)
old_down <- DustDistance(period_est_old[,1],tms,tem_old)
old_full <- DustDistance(period_est_old_FULL[,1],tms_FULL,tem_old)

## obtain true distances and dust from sesar table
sesar <- cbind(distance[1:N],extcr[1:N])


## plot point in black if period is correct, red if period is wrong
is_correct <- function(true,est) abs((true - est)/true) < .01
y1col <- (!is_correct(periods,period_est_new[,1])) + 1
y2col <- (!is_correct(periods,period_est_new_FULL[,1])) + 1
y3col <- (!is_correct(periods,period_est_old[,1])) + 1
y4col <- (!is_correct(periods,period_est_old_FULL[,1])) + 1

print("downsampled fraction of periods correct with new templates:")
mean(is_correct(periods,period_est_new[,1]))
print("downsampled fraction of periods correct with old templates:")
mean(is_correct(periods,period_est_old[,1]))

print("FULL fraction of periods correct with new templates:")
mean(is_correct(periods,period_est_new_FULL[,1]))
print("FULL fraction of periods correct with old templates:")
mean(is_correct(periods,period_est_old_FULL[,1]))




####### distance estimates
pdf("figs/1-results-distance.pdf",height=9,width=9)
ScatterMatrix(sesar[,1],new_down[,1],new_full[,1],
              old_down[,1],old_full[,1],
              xlab="Sesar Distance",
              r1lab=expression(bold("Downsampled Light Curve")),
              r2lab=expression(bold("Well Sampled Light Curve")),
              c1lab="Dependent Templates",
              c2lab="Fixed Templates",lim=c(5,120),
              y1col=y1col,y2col=y2col,y3col=y3col,y4col=y4col)
dev.off()

####### dust estimates
lim <- range(c(sesar[,2],new_down[,2],new_full[,2],old_down[,2],old_full[,2]))
pdf("figs/1-results-dust.pdf",height=9,width=9)
ScatterMatrix(sesar[,2],new_down[,2],new_full[,2],
              old_down[,2],old_full[,2],
              xlab="Schlegel r-band Extinction",
              r1lab=expression(bold("Downsampled Light Curve")),
              r2lab=expression(bold("Well Sampled Light Curve")),
              c1lab="Dependent Templates",
              c2lab="Fixed Templates",lim=lim,
              y1col=y1col,y2col=y2col,y3col=y3col,y4col=y4col)
dev.off()


kpc_to_mu <- function(kpc) 5*(log10(kpc*1000)-1)
mu_to_kpc <- function(mu) (10^(mu/5 + 1)) / 1000




plot(kpc_to_mu(new_down[,1]) - kpc_to_mu(sesar[,1]),new_down[,2]-sesar[,2])



y <- kpc_to_mu(new_full[,1]) - kpc_to_mu(sesar[,1])
x <- new_full[,2]-sesar[,2]
lm.fit <- lm(y~x)
plot(x,y)
abline(lm.fit$coeff)
lm.fit$coeff
abline(v=0)
abline(h=0)

pdf("extc_mu.pdf")
par(mar=c(5,5,1,1))
plot(x,y,
     ylab="mu Estimate (RR Model) - mu True (sesar)",xlab="Extinction r est (RR Model) - True Extinction r (Schlegel)",
     main="Well Sampled SDSS Stripe 82 RRL")
abline(lm.fit$coeff)
abline(v=0)
abline(h=0)
dev.off()




pred_mu <- kpc_to_mu(new_full[,1]) - predict(lm.fit)
pdf("distance_estimates_corrected.pdf",width=12,height=6)
par(mfcol=c(1,2),mar=c(5,5,1,1))
plot(sesar[,1],new_full[,1],xlab="Sesar Distance",ylab="Model Estimate")
plot(sesar[,1],mu_to_kpc(pred_mu),xlab="Sesar Distance",ylab="Model Estimate (Corrected)")
abline(a=0,b=1.01)
abline(a=0,b=0.99)
dev.off()



hist((mu_to_kpc(pred_mu) - sesar[,1])/sesar[,1])


par(mfcol=c(1,2))
hist(mu_to_kpc(pred_mu)/sesar[,1])
hist(new_full[,1]/sesar[,1])


## mui <- kpc_to_mu(new_full[,1])
## muis <- kpc_to_mu(distance[1:N])
## xi <- log10(periods + .2)
## xi2 <- xi^2


## lm.fit <- lm(muis - mui ~ xi + xi2)
## summary(mui-muis)
## hist(mui-muis)
## sd(mui-muis)
## summary(lm.fit$residuals)
## hist(lm.fit$residuals)
## sd(lm.fit$residuals)


## par(mfcol=c(1,2))
## plot(sesar[,1],new_full[,1],main="Original")
## abline(a=0,b=1)
## abline(a=0,b=1.05)
## abline(a=0,b=.95)
## lm.fit <- lm(muis - mui ~ 1)
## plot(mu_to_kpc(muis),mu_to_kpc(mui+predict(lm.fit)),main="Intercept Only")
## abline(a=0,b=1)
## abline(a=0,b=1.05)
## abline(a=0,b=.95)
