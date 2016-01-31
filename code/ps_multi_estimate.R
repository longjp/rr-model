rm(list=ls())
set.seed(1234)

library(parallel)

##source("func_multi_sine.R")
source("func_multi_saw.R")
source("func.R")

load("ps_multi_sim.RData")


period_min <- 0.2
period_max <- 1
mc.cores <- 7
NN <- 5
newton_pest <- function(ii){
    print(ii)
    lc <- tms[[ii]]
    omegas <- get_freqs(period_min,period_max,.1/diff(range(lc[,1])))
    omegas <- omegas / (2*pi)
    rss <- SawRss(m=lc$m,t=lc$t,ampj=lc$ampj,phij=lc$phij,omegas=omegas,NN=NN)
    return(1/omegas[which.min(rss)])
}
p_est_new <- unlist(mclapply(1:length(tms),newton_pest,mc.cores=mc.cores))


## K <- 1
## sine_pest <- function(ii){
##     print(ii)
##     lc <- tms[[ii]]
##     omegas <- get_freqs(period_min,period_max,.1/diff(range(lc[,1])))
##     lc <- to_list_lc(lc)
##     rss <- matrix(0,nrow=length(omegas),ncol=length(lc))
##     for(jj in 1:length(lc)){
##         rss[,jj] <- vapply(omegas,compute_rss,c(0),K,lc[[jj]])
##     }
##     rss <- rowSums(rss)
##     return((2*pi)/omegas[which.min(rss)])
## }
## p_est_sine <- unlist(mclapply(1:length(tms),sine_pest,mc.cores=mc.cores))


save(tms,periods,p_est_new,file="ps_multi_estimate.RData")
