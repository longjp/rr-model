rm(list=ls())
load("tms_params.RData")
load("make_template.RData")
load("estimate_params.RData")
source("rrab_fit.R")
source('func.R')

unlink("figs",recursive=TRUE)
dir.create("figs")

## report accuracies
N <- nrow(period_est)
print("accuracies, top 5:")
print(paste("1%:",mean(within_x(period_est,param$period[1:N],0.01))))
print(paste("0.1%:",mean(within_x(period_est,param$period[1:N],0.001))))
print(paste("0.01%:",mean(within_x(period_est,param$period[1:N],0.0001))))
print("")

print("accuracies, top period:")
period_est <- period_est[,1] ## just use best fit period
print(paste("1%:",mean(abs((param$period[1:N] - period_est)/param$period[1:N]) < 0.01)))
print(paste("0.1%:",mean(abs((param$period[1:N] - period_est)/param$period[1:N]) < 0.001)))
print(paste("0.01%:",mean(abs((param$period[1:N] - period_est)/param$period[1:N]) < 0.0001)))


## create template functions
temp_time <- seq(0,1,length.out=ncol(tem$templates))
tem$template_funcs <- list()
for(jj in 1:nrow(tem$templates)){
    tem$template_funcs[[jj]] <- approxfun(temp_time,tem$templates[jj,])
}
tem$templated_funcs <- list()
for(jj in 1:nrow(tem$templatesd)){
    tem$templated_funcs[[jj]] <- approxfun(temp_time,tem$templatesd[jj,])
}

## for a given omega, find coefficients
ComputeCoeffs <- function(tm,omega,tem,NN=10){
    dat <- AugmentData(tm,tem$dust,tem$betas)
    m <- dat[[1]]$mag
    dust <- dat[[1]]$dust
    t <- dat[[1]]$time
    nb <- dat[[2]]
    coeffs <- c(0,0,0,runif(1))
    while(coeffs[3]==0){
        for(jj in 1:NN){
            coeffs <- NewtonUpdate(coeffs[4],omega,m,t,dust,nb,tem$template_funcs,tem$templated_funcs)
        }
    }
    return(coeffs)
}



## plot all bands with best fit parameters, store in figs
for(ii in 1:N){
    tm <- tms[[ii]]
    omega <- 1/period_est[ii]
    coeffs <- ComputeCoeffs(tm,omega,tem)
    pdf(paste0("figs/",ii,".pdf"),height=12,width=8)
    dat <- SplitData(tm)
    par(mar=c(3,4,2,1),mfcol=c(5,1))
    for(jj in 1:length(tem$dust)){
        pred <- (coeffs[1] + tem$betas[jj] + coeffs[2]*tem$dust[jj]
            + coeffs[3]*tem$template_funcs[[jj]]((temp_time + coeffs[4]) %% 1))
        plot(temp_time/omega,pred,type='l',xlab="Phase",ylab="Mag")
        points((dat[[jj]]$time %% (1/omega)),dat[[jj]]$mag)
    }
    dev.off()
}
