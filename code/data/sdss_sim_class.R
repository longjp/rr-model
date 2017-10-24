## downsample sdss stripe 82 rr lyrae ab to Nobs total observations
rm(list=ls())
set.seed(1234)
source('../common/funcs.R')

####### get RRL light curves from table 1 provided by sesar
ToLC <- function(x){
    band <- rep(c("u","g","r","i","z"),each=nrow(x))
    names(x) <- rep(c("time","mag","error"),5)
    x <- rbind(x[,1:3],x[,4:6],x[,7:9],x[,10:12],x[,13:15])
    x <- data.frame(time=x[,1],band=band,mag=x[,2],error=x[,3])
    x[,2] <- as.character(x[,2])
    return(x[x[,3] > -90,])
}

fs <- list.files("raw/table1",full=TRUE)
fs <- fs[!grepl("README",fs)]

lcsRR <- vector("list",length(fs))
for(ii in 1:length(lcsRR)){
    x <- read.table(fs[ii])
    lcsRR[[ii]] <- ToLC(x[,3:ncol(x)])
}

lcsRR_fnames <- sub(".dat","",sub("raw/table1/","",fs))
names(lcsRR) <- lcsRR_fnames
lcsRR <- lcsRR[order(names(lcsRR))]


## load rrl and get distances
rrlyrae <- read.table("raw/apj326724t2_mrt.txt",skip=42)
temp <- read.table("raw/apj326724t3_mrt.txt",skip=30) ## get distances for rrlyrae, in different file
identical(rrlyrae$V1,temp$V1) ## ids same
rrlyrae$d <- temp$V5 ## put distances in rrlyrae data frame
rrlyrae$ra <- temp$V2
rrlyrae$dec <- temp$V3
rrlyrae$extcr <- temp$V4
rrlyrae[,1] <- as.character(rrlyrae[,1])
rrlyrae <- rrlyrae[order(rrlyrae[,1]),]
identical(rrlyrae[,1],names(lcsRR))
names(rrlyrae)[1] <- "ID"

## use only "ab" rrlyrae
rrab <- rrlyrae[,2] == "ab"
rrlyrae <- rrlyrae[rrab,]  ## only use rr lyrae ab
lcsRR <- lcsRR[rrab]


periodsRR <- rrlyrae[,3]
raRR <- rrlyrae$ra
decRR <- rrlyrae$dec
distanceRR <- rrlyrae$d
extcrRR <- rrlyrae$extcr

## load light curves
folder <- "raw/AllLCs/"
fnames <- list.files(folder)
fnames <- fnames[fnames!="LC_reorganize.tcl.dat"] ## get rid of non lc file

## only use well sampled light curves
nr <- vapply(fnames,function(x){nrow(read.table(paste(folder,x,sep="/")))},c(0))
fnames <- fnames[nr > 75]

## get non--rlyrae
Nnot <- 1000 ## number of non--rrlyrae to select
fnot <- fnames[!(fnames %in% paste0("LC_",lcsRR_fnames,".dat"))]
fnot <- sample(fnot,Nnot)
cat <- read.table("raw/stripe82candidateVar_v1.1.dat",header=TRUE)
ids <- gsub(".dat","",gsub("LC_","",fnot))
mean(ids %in% cat$ID) ## should = 1 (ie all ids in catalog)
raNOTRR <- cat$ra[order(match(cat$ID,ids))][1:length(ids)] # beautiful 1-liner
decNOTRR <- cat$dec[order(match(cat$ID,ids))][1:length(ids)] # beautiful 1-liner


## load dust
nf <- sum(grepl(".txt",list.files("extc/radec_extc")))
extcr <- vector("list",nf)
for(ii in 1:nf){
    extcr[[ii]] <- read.csv(paste0("extc/radec_extc/",ii-1,".txt"))
}
extcr <- do.call(rbind,extcr)
extcr[,1] <- as.character(cat$ID[1:nrow(extcr)])
names(extcr)[1] <- "ID"
extcr <- extcr[,c("ID","E_B_V_SFD")]
R_r <- 2.751 ## from dust website
extcr[,2] <- extcr[,2]*R_r
colnames(extcr)[2] <- "extcr"
extcrNOTRR <- extcr$extcr[order(match(extcr$ID,ids))][1:length(ids)] # beautiful 1-liner



## check that rrl ebv match
names(extcr)[2] <- "extcr_web"
out <- merge(extcr,rrlyrae[,c("ID","extcr")])
plot(out$extcr_web,out$extcr)
abline(a=0,b=1)



## QUESTION: HOW TO MERGE INTO ebv values

## compile attributes for each light curve
periods <- c(periodsRR,rep(0,Nnot)) ## if not rrl, period=0
distance <- c(distanceRR,rep(0,Nnot)) # if not rrl, distance=0
extcr <- c(extcrRR,extcrNOTRR)
ra <- c(raRR,raNOTRR)
dec <- c(decRR,decNOTRR)
cl <- c(rep("rr",length(lcsRR)),rep("not",Nnot))

## find ra, dec for all fnames light curves

## load light curves
lcsNOT <- vector("list",length(fnot))
for(ii in 1:length(fnot)){
    lcsNOT[[ii]] <- read.table(paste(folder,fnot[ii],sep="/"),stringsAsFactors=FALSE)
}
for(ii in 1:length(lcsNOT)) names(lcsNOT[[ii]]) <- c("time","band","mag","error")
names(lcsNOT) <- fnot



## merge nonRRL and RRL
lcs <- c(lcsRR,lcsNOT)


## save original data
tms_FULL <- lapply(lcs,LCtoTM)
names(tms_FULL) <- names(lcs)

## downsampled to total of Nobs observations
Nobs <- 20
for(ii in 1:length(lcs)){
    lc <- lcs[[ii]]
    lc <- lc[lc[,2]!="u",]
    ix <- sample(1:nrow(lc),Nobs)
    lcs[[ii]] <- lc[ix,]
}
tms <- lapply(lcs,LCtoTM)
names(tms) <- names(lcs)

## output lightcurves and periods
save(tms,tms_FULL,periods,cl,ra,dec,distance,extcr,file="clean/sdss_sim_class.RData")
