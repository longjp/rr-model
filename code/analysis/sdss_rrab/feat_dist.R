rm(list=ls())

load("../../fit_template/template_sdss.RData")
load("../../data/clean/sdss_rrab.RData")
source("../../fit_template/fit_template.R")
source("../../common/funcs.R")
source("../funcs.R")


lcs <- lapply(tms,TMtoLC)

## find the residual sum of squares from model fit
N <- length(lcs)
rss.n <- rep(0,N)
coeffs <- matrix(0,nrow=length(lcs),ncol=4)
for(ii in 1:N){
    print(ii)
    lc <- lcs[[ii]]
    omega <- 1/periods[ii]
    coes <- ComputeCoeffs(lc,omega,tem_sdss,NN=20)
    coeffs[ii,] <- coes
    pred <- PredictTimeBand(lc$time,lc$b,omega,coes,tem_sdss)
    rss.n[ii] <- median(abs(pred - lc$mag))
}


###### TODO:
## 1. estimate distances old way, compare to sesar, should see 10% scatter
## 2. use sesar 2010 and estimated amplitudes to improve distance estimate
## 3. how to present results, is this somehow accounting for period?


coeff <- cbind(coeffs,periods,rss.n)
colnames(coeff) <- c("mu","ebv","amp_g (p_to_p)","phi","period","MAR")
##png("stripe82.png",height=600,width=600)
to_use <- (coeff[,"amp_g (p_to_p)"] < 2) & (coeff[,"ebv"] < .75)  & (coeff[,"ebv"] > -.25) & (coeff[,"MAR"] < .15)## get rid of extreme values
cols <- c("red")
pchs <- c(2)
columns <- c("amp_g (p_to_p)","period","MAR")
par(cex.axis=1.5)
pairs(coeff[to_use,columns],col=cols,pch=pchs,
      main="Sloan Stripe 82")
##dev.off()

par(mar=c(4.5,4.5,1,1))
plot(coeff[to_use,5],coeff[to_use,3],col='red',pch=2,xlim=c(.2,1),ylim=c(0,2),xlab="Period",ylab="amp_g (p_to_p)",
     cex.lab=1.5)


head(coeff)
plot(coeff[,2],coeff[,5],xlab="Dust",ylab="period")
lm.fit <- lm(coeff[,5]~coeff[,2])
abline(lm.fit)
plot(coeff[,3],coeff[,5],xlab="Amp",ylab="period")


## load dust and distances for these stars
## plot distance residuals (sesar distance - our distance) versus period, amplitude


## plot dust residual (sesar dust - out dust) versus period

