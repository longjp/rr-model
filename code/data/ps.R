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

### checking if TMtoLC and LCtoTM are inverses, issues with ordering
### of measurements. verify that these are inverses
## lcs2 <- lapply(tms,TMtoLC)
## identical(lcs,lcs2)
## lcs[[1]]
## lcs2[[1]]
## lcs[[1]][,2]==lcs2[[1]][,2]
## lcs[[1]][,2]
## lcs2[[1]][,2]

## output lightcurves and periods
save(tms,periods,file="clean/ps.RData")
