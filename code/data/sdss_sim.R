## downsample sdss stripe 82 rr lyrae ab to 50 total observations
rm(list=ls())
source('funcs.R')

rrlyrae <- read.table("raw/apj326724t2_mrt.txt",skip=42)
rrlyrae <- rrlyrae[rrlyrae[,2] == "ab",]

folder <- "raw/AllLCs/"
fnames <- list.files(folder)

## order fnames and rrlyrae
rrlyrae[,1] <- paste("LC_",rrlyrae[,1],".dat",sep="")
fnames <- fnames[fnames %in% rrlyrae[,1]]
rrlyrae <- rrlyrae[rrlyrae[,1] %in% fnames,]
rrlyrae <- rrlyrae[order(rrlyrae[,1]),]
fnames <- fnames[order(fnames)]
periods <- rrlyrae[,3]

## load light curves
lcs <- list()
for(ii in 1:length(fnames)){
    lcs[[ii]] <- read.table(paste(folder,fnames[ii],sep="/"))
}
for(ii in 1:length(lcs)) names(lcs[[ii]]) <- c("time","band","mag","sigma")

## downsampled to total of Nobs observations
Nobs <- 50
for(ii in 1:length(lcs)){
    ix <- sample(1:nrow(lcs[[ii]]),Nobs)
    lcs[[ii]] <- lcs[[ii]][ix,]
}
tms <- lapply(lcs,LCtoTMS)
names(tms) <- rrlyrae[,1]

## output lightcurves and periods
save(tms,periods,file="clean/sdss_sim.RData")

