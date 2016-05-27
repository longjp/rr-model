rm(list=ls())

## load necessary libraries
library('parallel')
library('multiband')
load("../../fit_template/template.RData")
source("../../fit_template/template.R")
source("../../common/funcs.R")
source("../funcs.R")
source("../params.R")

## data source
load("../../data/clean/sdss_sim.RData")
load("results.RData")

period_est <- period_est[,1] ## just use best fit period


tot.dev <- rep(0,N)
for(ii in 1:N){
    tm <- tms[[ii]]
    omega <- 1/period_est[ii]
    lc <- TMtoLC(tm)
    coeffs <- ComputeCoeffs(lc,omega,tem)
    dev <- rep(0,length(tm))
    for(jj in 1:length(tm)){
        pred <- (coeffs[1] + tem$betas[jj] + coeffs[2]*tem$dust[jj]
            + coeffs[3]*tem$template_funcs[[jj]]((tm[[jj]][,1]*omega + coeffs[4]) %% 1))
        dev[jj] <- sum(abs((pred - tm[[jj]][,2])))
    }
    tot.dev[ii] <- sum(dev)
}
