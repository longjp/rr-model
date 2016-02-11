## construct .RData file for ps_multi_estimate.R to analyze
## output is tms (from panstarrs), periods (truth), cc = fraction down
rm(list=ls())

load("sdss_eda.RData")


rrlyrae <- read.table("apj326724t2_mrt.txt",skip=42)
rrlyrae <- rrlyrae[rrlyrae[,2] == "ab",]

folder <- "rrlyrae"
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
for(ii in 1:length(tms)) names(tms[[ii]]) <- c("t","b","m","sigma")


for(ii in 1:length(tms)){
    band <- as.character(tms[[ii]]$b)
    tms[[ii]] <- data.frame(t=tms[[ii]][,1],
                            m=tms[[ii]][,3]-betas[band],
                            sig=1,ampj=amps[band],
                            phij=phis[band])
}


names(tms) <- rrlyrae[,1]
periods <- rrlyrae[,3]

save(tms,periods,cc,file="ps_multi.RData")
