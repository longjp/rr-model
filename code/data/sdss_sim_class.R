## downsample sdss stripe 82 rr lyrae ab to Nobs total observations
rm(list=ls())
set.seed(1234)
source('../common/funcs.R')

## load rrl and get distances
rrlyrae <- read.table("raw/apj326724t2_mrt.txt",skip=42)
temp <- read.table("raw/apj326724t3_mrt.txt",skip=30) ## get distances for rrlyrae, in different file
identical(rrlyrae$V1,temp$V1) ## ids same
rrlyrae$d <- temp$V5 ## put distances in rrlyrae data frame

## only use rr lyrae ab
rrlyrae <- rrlyrae[rrlyrae[,2] == "ab",]

## load light curves
folder <- "raw/AllLCs/"
fnames <- list.files(folder)
fnames <- fnames[fnames!="LC_reorganize.tcl.dat"] ## get rid of non lc file

## only use well sampled light curves
nr <- vapply(fnames,function(x){nrow(read.table(paste(folder,x,sep="/")))},c(0))
fnames <- fnames[nr > 75]

## order fnames and rrlyrae
Nnot <- 1000 ## number of non--rrlyrae to select
rrlyrae[,1] <- paste0("LC_",rrlyrae[,1],".dat")
rrlyrae <- rrlyrae[order(rrlyrae[,1]),]
fnot <- sample(fnames[!(fnames %in% rrlyrae[,1])],Nnot)
fnames <- fnames[fnames %in% rrlyrae[,1]]
fnames <- fnames[order(fnames)]
rrlyrae <- rrlyrae[rrlyrae[,1] %in% fnames,]
fnames <- c(fnames,fnot)
periods <- c(rrlyrae[,3],rep(0,Nnot))
distance <- c(rrlyrae$d,rep(0,Nnot))
cl <- c(rep("rr",nrow(rrlyrae)),rep("not",Nnot))

## find ra, dec for all fnames light curves
cat <- read.table("raw/stripe82candidateVar_v1.1.dat",header=TRUE)
ids <- gsub(".dat","",gsub("LC_","",fnames))
mean(ids %in% cat$ID) ## should = 1 (ie all ids in catalog)
ra <- cat$ra[order(match(cat$ID,ids))][1:length(ids)] # beautiful 1-liner
dec <- cat$dec[order(match(cat$ID,ids))][1:length(ids)] # beautiful 1-liner

## load light curves
lcs <- vector("list",length(fnames))
for(ii in 1:length(fnames)){
    lcs[[ii]] <- read.table(paste(folder,fnames[ii],sep="/"))
}
for(ii in 1:length(lcs)) names(lcs[[ii]]) <- c("time","band","mag","sigma")


## save original data
tms_FULL <- lapply(lcs,LCtoTM)
names(tms_FULL) <- fnames

## downsampled to total of Nobs observations
Nobs <- 20
for(ii in 1:length(lcs)){
    lc <- lcs[[ii]]
    lc <- lc[lc[,2]!="u",]
    ix <- sample(1:nrow(lc),Nobs)
    lcs[[ii]] <- lc[ix,]
}
tms <- lapply(lcs,LCtoTM)
names(tms) <- fnames

## output lightcurves and periods
save(tms,tms_FULL,periods,cl,ra,dec,distance,file="clean/sdss_sim_class.RData")
