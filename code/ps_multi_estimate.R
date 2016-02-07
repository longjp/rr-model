rm(list=ls())
set.seed(1234)

library(parallel)
library(multiband)

source("func_multi_saw.R")
source("func.R")

load("ps_multi.RData")

period_min <- 0.2
period_max <- 1
mc.cores <- 7

print("running Newton NN=1")
NN <- 1
newton_pest <- function(ii){
    print(ii)
    lc <- tms[[ii]]
    omegas <- get_freqs(period_min,period_max,.1/diff(range(lc[,1])))
    omegas <- omegas / (2*pi)
    rss <- SawRss(m=lc$m,t=lc$t,ampj=lc$ampj,phij=lc$phij,omegas=omegas,NN=NN,cc=cc)
    return(1/omegas[sort_local_min(1:length(rss),rss)[1:5]])
}
p_est_new1 <- mclapply(1:length(tms),newton_pest,mc.cores=mc.cores)
p_est_new1 <- matrix(unlist(p_est_new1),ncol=5,byrow=TRUE)



print("running Newton NN=5")
NN <- 5
newton_pest <- function(ii){
    print(ii)
    lc <- tms[[ii]]
    omegas <- get_freqs(period_min,period_max,.1/diff(range(lc[,1])))
    omegas <- omegas / (2*pi)
    rss <- SawRss(m=lc$m,t=lc$t,ampj=lc$ampj,phij=lc$phij,omegas=omegas,NN=NN,cc=cc)
    return(1/omegas[sort_local_min(1:length(rss),rss)[1:5]])
}
p_est_new5 <- unlist(mclapply(1:length(tms),newton_pest,mc.cores=mc.cores))
p_est_new5 <- matrix(unlist(p_est_new5),ncol=5,byrow=TRUE)


print("running Newton NN=10")
NN <- 10
newton_pest <- function(ii){
    print(ii)
    lc <- tms[[ii]]
    omegas <- get_freqs(period_min,period_max,.1/diff(range(lc[,1])))
    omegas <- omegas / (2*pi)
    rss <- SawRss(m=lc$m,t=lc$t,ampj=lc$ampj,phij=lc$phij,omegas=omegas,NN=NN,cc=cc)
    return(1/omegas[sort_local_min(1:length(rss),rss)[1:5]])
}
p_est_new10 <- unlist(mclapply(1:length(tms),newton_pest,mc.cores=mc.cores))
p_est_new10 <- matrix(unlist(p_est_new10),ncol=5,byrow=TRUE)


print("running sine")
sine_pest <- function(ii){
    print(ii)
    omegas <- get_freqs(period_min,period_max,.1/diff(range(tms[[ii]][,1])))
    lc <- split(tms[[ii]],tms[[ii]]$ampj)
    out <- pgls(lc,periods = (2*pi)/omegas,BCD_flag=FALSE)
    ix <- sort_local_min(1:length(out$period_seq_all),out$rss_ls)
    return(out$period_seq_all[ix[1:5]])
}
p_est_sine <- mclapply(1:length(tms),sine_pest,mc.cores=mc.cores)
p_est_sine <- matrix(unlist(p_est_sine),ncol=5,byrow=TRUE)


save(tms,periods,p_est_new1,p_est_new5,p_est_new10,p_est_sine,file="ps_multi_estimate.RData")
