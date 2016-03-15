## code simulate RR Lyrae from model with sdss cadence
rm(list=ls())
load("make_template.RData")
load("make_tms.RData")

bands <- names(amps)
t <- (0:ncol(templates))/ncol(templates)

## simulate parameters all at once
N <- 500 # number of light curves to simulate
alpha <- rnorm(N)
d <- rnorm(N)
a <- rep(0.3,N)
period <- runif(N,min=0.4,max=0.8)
phase <- runif(N,min=0,max=1)

## store output
tms_sim <- list()

## simulate light curves
for(ii in 1:N){
    ix <- sample(1:length(tms),1)
    tms_sim[[ii]] <- tms[[ix]]
    tms_sim[[ii]][,3] <- 0
    t_temp <- (tms_sim[[ii]][,1] %% period[ii]) / period[ii]
    t_temp <- (t_temp + phase[ii]) %% 1
    for(jj in 1:length(bands)){
        to_use <- tms_sim[[ii]][,2] == bands[jj]
        m <- approx(t,c(templates[jj,ncol(templates)],templates[jj,]),xout=t_temp[to_use])$y
        m <- betas[jj] + alpha[ii] + d[ii]*dust[jj] + a[ii]*amps[jj]*m + 3*rnorm(length(m),mean=0,sd=tms_sim[[ii]][,4])
        tms_sim[[ii]][to_use,3] <- m
    }
}


## ii <- 200
## plot((tms_sim[[ii]][,1] %% period[ii])/period[ii],tms_sim[[ii]][,3],col=as.numeric(tms_sim[[ii]][,2]))

save(tms_sim,alpha,a,period,d,phase,file="sim.RData")
