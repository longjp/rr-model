library('parallel')
source("rrab_fit.R")
load("tms_params.RData")
load("make_template.RData")

temp_time <- seq(0,1,length.out=ncol(tem$templates))
tem$template_funcs <- list()
for(jj in 1:nrow(tem$templates)){
    tem$template_funcs[[jj]] <- approxfun(temp_time,tem$templates[jj,])
}
tem$templated_funcs <- list()
for(jj in 1:nrow(tem$templatesd)){
    tem$templated_funcs[[jj]] <- approxfun(temp_time,tem$templatesd[jj,])
}

ComputePeriod <- function(ii){
    print(ii)
    rss <- ComputeRSS(tms[[ii]],omegas,tem,NN=NN)
    return(1/omegas[which.min(rss)])
}

## parameters for simulation
##N <- length(tms) ## number of light curves to run
N <- 3
NN <- 4
omegas <- get_freqs(0.2,1)
mc.cores <- 2

## estimate periods
period_est <- mclapply(1:N,ComputePeriod,mc.cores=mc.cores)
period_est <- as.numeric(period_est)

save(period_est,file="estimate_params.RData")
