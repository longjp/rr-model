rm(list=ls())
source('func.R')
library('multiband')
library('parallel')
load("tms_params.RData")

ComputePeriod <- function(ii){
    print(ii)
    tm <- tms[[ii]]
    bands <- unique(tm$band)
    lc <- list()
    for(jj in 1:length(bands)){
        if(sum(tm$band==bands[jj]) > 0.5){
            lc[[jj]] <- cbind(tm[tm$band==bands[jj],c(1,3)],1)
        }
    }
    out <- pgls(lc,periods=periods,BCD_flag=FALSE)
    rss <- out$rss_ls
    ps <- periods[sort_local_min(1:length(rss),rss)]
    ps <- SeparateBest(ps,0.0002,5)
    return(ps)
}


N <- 2
periods <- (2*pi) / get_freqs(0.2,1)
mc.cores <- 1

period_est_lomb <- mclapply(1:N,ComputePeriod,mc.cores=mc.cores)
period_est_lomb <- matrix(unlist(period_est_lomb),ncol=5,byrow=TRUE)

save(period_est_lomb,file="estimate_params_lomb.RData")
