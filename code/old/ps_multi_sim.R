## appears to be period evolution of lc shape for rrab
rm(list=ls())

set.seed(1234)

amp <- 0.3 ## typical amplitude for rr lyrae
cc <- 0.25 ## fraction lightcurve mag goes down


######## get rr lyrae periods
stars <- read.table("apj326724t2_mrt.txt",skip=42)
periods <- stars[,3]


######## USE TIMES AND ERRORS FROM PANSTARRS DATA TO GENERATE SAWTOOTHS
## read data
f <- list.files("PS1_sample_LCs",full.names=TRUE)
tms <- list()
for(ii in 1:length(f)){
    temp <- read.table(f[ii],header=TRUE)
    t <- temp[,3]
    e.sd <- temp[,1]
    band <- temp[,2]
    B <- length(unique(band))
    amp_mag <- (0:(B-1))/10 + 1
    names(amp_mag) <- unique(band)
    ampj <- amp_mag[band]
    phi_mag <- (0:(B-1))/10
    names(phi_mag) <- unique(band)
    phij <- phi_mag[band]
    tms[[ii]] <- data.frame(t=t,m=0,sig=e.sd,ampj=ampj,phij=phij)
}


gammaf <- function(t,cc=0.25){
    t <- t %% 1
    return(((-2/cc)*t + 1)*(t < cc) + ((2/(1-cc))*t - (1+cc)/(1-cc))*(t >= cc))
}

sinef <- function(t){
    t <- t %% 1
    return(sin(2*pi*t))
}


## simulate magnitudes, replace errors with 1
periods <- sample(periods,size=length(tms),replace=TRUE)
for(ii in 1:length(tms)){
    phi <- runif(1,min=0,max=1)
    tms[[ii]][,2] <- (amp*tms[[ii]]$ampj*gammaf((tms[[ii]]$t/periods[ii]) + tms[[ii]]$phij + phi,cc)
                      + rnorm(nrow(tms[[ii]]),mean=0,sd=tms[[ii]]$sig))
    tms[[ii]][,3] <- 1
}

######## 3. WRITE OUT RESULTS
save(tms,periods,cc,file="ps_multi.RData")
