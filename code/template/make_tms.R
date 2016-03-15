rm(list=ls())

f <- list.files("../rrlyrae/",full.names=TRUE)
tms <- list()
for(ii in 1:length(f)){
    tms[[ii]] <- read.table(f[ii])
    names(tms[[ii]]) <- c("time","band","mag","error")
}

save(tms,file="make_tms.RData")
