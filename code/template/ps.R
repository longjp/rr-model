## construct .RData file for ps_multi_estimate.R to analyze
## output is tms (from panstarrs), periods (truth), cc = fraction down
rm(list=ls())

##load("sdss_eda.RData")
period <- read.table("../period_comparisons.txt",header=TRUE)
period[,1] <- as.character(period[,1])
period <- period[period$RRtype=="ab",]

fs <- list.files("../PS1_sample_LCs")
tms <- list()
for(ii in 1:length(fs)){
    temp <- read.table(paste0("../PS1_sample_LCs/",fs[ii]),header=TRUE)
    band <- as.character(temp$filterid)
    band[band=="y"] = "z"
    band <- factor(band,levels=c("g","i","r","u","z"))
    tms[[ii]] <- data.frame(time=temp[,3],band=band,mag=temp[,4],error=1)
}
names(tms) <- fs


## use only sources for which we have BOTH lc and period
period <- period[period[,1] %in% names(tms),]
period <- period[order(period[,1]),]
tms <- tms[names(tms) %in% period[,1]]
tms <- tms[order(names(tms))]
period <- period[,2]

param <- list(period=period)

save(tms,param,file="tms_params.RData")




