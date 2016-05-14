## construct .RData file for ps_multi_estimate.R to analyze
## output is tms (from panstarrs), periods (truth), cc = fraction down
rm(list=ls())

rrlyrae <- read.table("apj326724t2_mrt.txt",skip=42)
rrlyrae <- rrlyrae[rrlyrae[,2] == "ab",]

folder <- "../rrlyrae"
fnames <- list.files(folder)

## order fnames and rrlyrae
rrlyrae[,1] <- paste("LC_",rrlyrae[,1],".dat",sep="")
fnames <- fnames[fnames %in% rrlyrae[,1]]
rrlyrae <- rrlyrae[rrlyrae[,1] %in% fnames,]
rrlyrae <- rrlyrae[order(rrlyrae[,1]),]
fnames <- fnames[order(fnames)]


## load light curves and put in tmss format
tms <- list()
for(ii in 1:length(fnames)){
    tms[[ii]] <- read.table(paste(folder,fnames[ii],sep="/"))
}
for(ii in 1:length(tms)) names(tms[[ii]]) <- c("time","band","mag","error")

## downsampled to total of Nobs observations
Nobs <- 40
for(ii in 1:length(tms)){
    ix <- sample(1:nrow(tms[[ii]]),Nobs)
    tms[[ii]] <- tms[[ii]][ix,]
}

## output lightcurves and periods
names(tms) <- rrlyrae[,1]
period <- rrlyrae[,3]
param <- list(period=period)
save(tms,param,file="tms_params.RData")

