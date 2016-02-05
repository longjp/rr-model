## construct .RData file for ps_multi_estimate.R to analyze
## output is tms (from panstarrs), periods (truth), cc = fraction down
rm(list=ls())

load("sdss_eda.RData")
periods <- read.table("period_comparisons.txt",header=TRUE)
periods[,1] <- as.character(periods[,1])

fs <- list.files("PS1_sample_LCs")
tms <- list()
for(ii in 1:length(fs)){
    temp <- read.table(paste0("PS1_sample_LCs/",fs[ii]),header=TRUE)
    band <- as.character(temp$filterid)
    band[band=="y"] = "z"
    tms[[ii]] <- data.frame(t=temp[,3],m=temp[,4]-betas[band],sig=1,ampj=amps[band],phij=phis[band])
}
names(tms) <- fs


## use only sources for which we have BOTH lc and period
periods <- periods[periods[,1] %in% names(tms),]
periods <- periods[order(periods[,1]),]
tms <- tms[names(tms) %in% periods[,1]]
tms <- tms[order(names(tms))]
periods <- periods[,2]


save(tms,periods,cc,file="ps_multi.RData")




