## appears to be period evolution of lc shape for rrab
rm(list=ls())

set.seed(1234)

amp <- 0.3 ## typical amplitude for rr lyrae

######## get rr lyrae periods
stars <- read.table("apj326724t2_mrt.txt",skip=42)
periods <- stars[,3]


######## USE TIMES AND ERRORS FROM PANSTARRS DATA TO GENERATE SAWTOOTHS
## read data
f <- list.files("PS1_sample_LCs",full.names=TRUE)
tms <- list()
for(ii in 1:length(f)){
    temp <- read.table(f[ii],header=TRUE)
    tms[[ii]] <- data.frame(temp[,3],0,temp[,1])
    names(tms[[ii]]) <- c("time","mag","error")
}
summary(vapply(tms,nrow,0))

gammaf <- function(t,cc=0.75){
    t <- t %% 1
    return(((-2/cc)*t + 1)*(t < cc) + ((2/(1-cc))*t - (1+cc)/(1-cc))*(t >= cc))
}

sinef <- function(t){
    t <- t %% 1
    return(sin(2*pi*t))
}


## simulate magnitudes
periods <- sample(periods,size=length(tms),replace=TRUE)
for(ii in 1:length(tms)){
    phi <- runif(1,min=0,max=1)
    tms[[ii]][,2] <- (amp*gammaf((tms[[ii]][,1]/periods[ii]) + phi)
                      + rnorm(nrow(tms[[ii]]),mean=0,sd=tms[[ii]][,3]))
    tms[[ii]][,3] <- 1
}

######## 3. WRITE OUT RESULTS
save(tms,periods,file="ps_sim.RData")
