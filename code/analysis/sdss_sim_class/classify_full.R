## develop a classifier on well sampled light curves
## which has very good performance
## derive features other than model outputs


rm(list=ls())
set.seed(1234)

## load necessary libraries
library('parallel')
library('multiband')
library('randomForest')
load("../../fit_template/template.RData")
source("../../fit_template/fit_template.R")
source("../../common/funcs.R")
source("../funcs.R")


## data source
load("../../data/clean/sdss_sim_class.RData")
load("results.RData")
source("../params.R")








ps <- period_est_FULL[,1] ## just use best fit period


coeffs <- matrix(0,nrow=length(ps),ncol=4)
for(ii in 1:nrow(coeffs)){
    omega <- 1/ps[ii]
    coeffs[ii,] <- ComputeCoeffs(TMtoLC(tms_FULL[[ii]]),omega,tem,NN=20,use.errors=FALSE)
}
colnames(coeffs) <- c("mu","ebv","p2p-gband","phase")



to_use <- cl=="rr"

par(mfcol=c(1,2))
plot(ps[to_use],coeffs[to_use,3],ylim=c(0,2))
plot(periods[to_use],ps[to_use],ylim=c(.2,1),xlim=c(.2,1))
abline(a=0,b=1)


## we get 90% of period estimates correct on full light curves, doesn't seem so good
e <- 20 / (60*60*24)
table(abs(ps[to_use]-periods[to_use]) < e)


## get RSS measures, fit RF, get classifier accuracy, example lcs which are misclassified
## possible reasons for missclassification
##   1. rrl with incorrect parameter estimates
##   2. lcs labeled as non-rrl are actually rrl (or look a lot like them)
##   3. need more features


### TODOs

## possible reasons: 1) error distribution of l.c.s with wrong periods has a few wild outliers
##                   2) l.c.s. are poorly sampled
##                   3) code note finding model global min
