rm(list=ls())
source("../common/funcs.R")

cat <- read.table("raw/period_comparisons.txt",header=TRUE)
fs <- cat[,1]
periods <- cat[,2]

## read in all light curves, pretend y band is z band
lcs <- vector("list",length(fs))
for(ii in 1:length(fs)){
    lc <- read.table(paste0("raw/PS1_sample_LCs/",fs[ii]),
                     header=TRUE,stringsAsFactors=FALSE)
    lc <- lc[,c(3,2,4,1)]
    names(lc) <- c("time","band","mag","error")
    lc$band[lc$band=="y"] <- "z"
    lcs[[ii]] <- lc
}

tms <- lapply(lcs,LCtoTM)

## output lightcurves and periods
save(tms,periods,file="clean/ps.RData")



