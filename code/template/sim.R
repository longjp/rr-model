## code simulate RR Lyrae from model with sdss cadence
rm(list=ls())
load("make_template.RData")

f <- list.files("../rrlyrae/",full.names=TRUE)
tms <- list()
for(ii in 1:length(f)){
    tms[[ii]] <- read.table(f[ii])
    names(tms[[ii]]) <- c("time","band","mag","error")
}


bands <- names(dust)
t <- seq(0,1,length.out=ncol(templates))

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
        m <- approx(t,templates[jj,],xout=t_temp[to_use])$y
        m <- betas[jj] + alpha[ii] + d[ii]*dust[jj] + a[ii]*m + 3*rnorm(length(m),mean=0,sd=tms_sim[[ii]][,4])
        tms_sim[[ii]][to_use,3] <- m
    }
}


## ii <- 200
## plot((tms_sim[[ii]][,1] %% period[ii])/period[ii],tms_sim[[ii]][,3],col=as.numeric(tms_sim[[ii]][,2]))

tms <- tms_sim
param <- list(period=period,alpha=alpha,a=a,d=d,phase=phase)
save(tms,param,file="sim.RData")
