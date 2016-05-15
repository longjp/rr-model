library('parallel')
source("rrab_fit.R")
source("func.R")
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
    rss <- FitTemplate(tms[[ii]],omegas,tem,NN=NN)
    ps <- 1/omegas[sort_local_min(1:length(rss),rss)]
    ps <- SeparateBest(ps,0.0002,5)
    return(ps)
}

## parameters for simulation
N <- length(tms)
NN <- 10
omegas <- get_freqs(0.2,1)
mc.cores <- 12

## estimate periods
period_est <- mclapply(1:N,ComputePeriod,mc.cores=mc.cores)
period_est <- matrix(unlist(period_est),ncol=5,byrow=TRUE)

save(period_est,file="estimate_params.RData")
