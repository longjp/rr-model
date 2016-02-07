rm(list=ls())
set.seed(1234)

library(parallel)
library(multiband)

##source("func_multi_sine.R")
source("func_multi_saw.R")
source("func.R")

load("ps_multi.RData")

period_min <- 0.2
period_max <- 1
mc.cores <- 7




## NN <- 5
## ii <- 22
## lc <- tms[[ii]]
## omegas <- get_freqs(period_min,period_max,.1/diff(range(lc[,1])))
## omegas <- omegas / (2*pi)
## rss <- SawRss(m=lc$m,t=lc$t,ampj=lc$ampj,phij=lc$phij,omegas=omegas,NN=NN,cc=cc)

## plot(omegas,rss)
## abline(v=omegas[which.min(rss)],col='red',lwd=2)
## abline(v=1/periods[ii])

## bad cases
## 8, 9, 15, 20, 21


print("running Newton NN=1")
NN <- 1
newton_pest <- function(ii){
    print(ii)
    lc <- tms[[ii]]
    omegas <- get_freqs(period_min,period_max,.1/diff(range(lc[,1])))
    omegas <- omegas / (2*pi)
    rss <- SawRss(m=lc$m,t=lc$t,ampj=lc$ampj,phij=lc$phij,omegas=omegas,NN=NN,cc=cc)
    return(1/omegas[which.min(rss)])
}
p_est_new1 <- unlist(mclapply(1:length(tms),newton_pest,mc.cores=mc.cores))

print("running Newton NN=5")
NN <- 5
newton_pest <- function(ii){
    print(ii)
    lc <- tms[[ii]]
    omegas <- get_freqs(period_min,period_max,.1/diff(range(lc[,1])))
    omegas <- omegas / (2*pi)
    rss <- SawRss(m=lc$m,t=lc$t,ampj=lc$ampj,phij=lc$phij,omegas=omegas,NN=NN,cc=cc)
    return(1/omegas[which.min(rss)])
}
p_est_new5 <- unlist(mclapply(1:length(tms),newton_pest,mc.cores=mc.cores))


print("running Newton NN=10")
NN <- 10
newton_pest <- function(ii){
    print(ii)
    lc <- tms[[ii]]
    omegas <- get_freqs(period_min,period_max,.1/diff(range(lc[,1])))
    omegas <- omegas / (2*pi)
    rss <- SawRss(m=lc$m,t=lc$t,ampj=lc$ampj,phij=lc$phij,omegas=omegas,NN=NN,cc=cc)
    return(1/omegas[which.min(rss)])
}
p_est_new10 <- unlist(mclapply(1:length(tms),newton_pest,mc.cores=mc.cores))


print("running sine")
sine_pest <- function(ii){
    print(ii)
    omegas <- get_freqs(period_min,period_max,.1/diff(range(tms[[ii]][,1])))
    lc <- split(tms[[ii]],tms[[ii]]$ampj)
    out <- pgls(lc,periods = (2*pi)/omegas,BCD_flag=FALSE)
    return(out$period_seq_all[which.min(out$rss_ls)])
}
p_est_sine <- unlist(mclapply(1:length(tms),sine_pest,mc.cores=mc.cores))


save(tms,periods,p_est_new1,p_est_new5,p_est_new10,p_est_sine,file="ps_multi_estimate.RData")
