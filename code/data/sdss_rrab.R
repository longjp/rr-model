## select RR ab for template creation, template creation selects well sampled
rm(list=ls())
source('../common/funcs.R')

rrlyrae <- read.table("raw/apj326724t2_mrt.txt",skip=42)
temp <- read.table("raw/apj326724t3_mrt.txt",skip=30) ## get distances for rrlyrae, in different file
identical(rrlyrae$V1,temp$V1) ## ids same
rrlyrae$d <- temp$V5 ## put distances in rrlyrae data frame
rrlyrae$ra <- temp$V2
rrlyrae$dec <- temp$V3
rrlyrae$extcr <- temp$V4
rrlyrae <- rrlyrae[rrlyrae[,2] == "ab",]

folder <- "raw/AllLCs"
fnames <- list.files(folder)

## order fnames and rrlyrae, get periods, other meta data
rrlyrae[,1] <- paste("LC_",rrlyrae[,1],".dat",sep="")
fnames <- fnames[fnames %in% rrlyrae[,1]]
rrlyrae <- rrlyrae[rrlyrae[,1] %in% fnames,]
rrlyrae <- rrlyrae[order(rrlyrae[,1]),]
fnames <- fnames[order(fnames)]
periods <- rrlyrae[,3]
ra <- rrlyrae$ra
dec <- rrlyrae$dec
distance <- rrlyrae$d
extcr <- rrlyrae$extcr

## load light curves and put in tms format
lcs <- vector("list",length(fnames))
for(ii in 1:length(fnames)){
    lcs[[ii]] <- read.table(paste(folder,fnames[ii],sep="/"),stringsAsFactors=FALSE)
}
for(ii in 1:length(lcs)) names(lcs[[ii]]) <- c("time","band","mag","error")
tms <- lapply(lcs,LCtoTM)
names(tms) <- rrlyrae[,1]

## save
save(tms,periods,ra,dec,distance,extcr,file="clean/sdss_rrab.RData")
