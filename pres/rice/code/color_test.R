## make feature plots showing rrl versus nonrrl separation
rm(list=ls())
library(scales)

fold <- "../../../code/data/raw/AllLCs/"
cat <- read.table("stripe82candidateVar_v1.1.dat",header=TRUE)
cat[,1] <- paste0("LC_",cat[,1],".dat")
fs <- list.files(fold)
fs <- fs[fs!="LC_reorganize.tcl.dat"] ## get rid of non lc file

cat <- cat[cat[,1] %in% fs,]
fs <- fs[fs %in% cat[,1]]
cat <- cat[order(cat[,1]),]
fs <- fs[order(fs)]

## read in rrlyrae
rrlyrae <- read.table("apj326724t2_mrt.txt",skip=42)
rrlyrae <- rrlyrae[paste0("LC_",rrlyrae[,1],".dat") %in% fs,]
rrlyrae[,1] <- paste0("LC_",rrlyrae[,1],".dat")
## period, ug color, amplitude
rrlyrae <- rrlyrae[,c(1,2,3,8)]
colnames(rrlyrae) <- c("ID","subtype","p","gamp")

## ComputeColor <- function(lc,b1,b2){
##     return(median(lc[lc[,2]==b1,3]) - median(lc[lc[,2]==b2,3]))
## }

## lcs <- vector("list",length(fs))
## for(ii in 1:length(fs)){
##     lcs[[ii]] <- read.table(paste0(fold,fs[ii]))
## }

## ## compute colors
## b1 <- "u"
## b2 <- "g"
## ug <- vapply(lcs,ComputeColor,c(0),b1=b1,b2=b2)
## b1 <- "g"
## b2 <- "r"
## gr <- vapply(lcs,ComputeColor,c(0),b1=b1,b2=b2)
## b1 <- "r"
## b2 <- "i"
## ri <- vapply(lcs,ComputeColor,c(0),b1=b1,b2=b2)
## b1 <- "i"
## b2 <- "z"
## iz <- vapply(lcs,ComputeColor,c(0),b1=b1,b2=b2)

## merge catalog with rrlyrae, use sesar feature estimates
cat <- merge(cat,rrlyrae,all=TRUE)
cat$P[!is.na(cat$p)] <- cat$p[!is.na(cat$p)]
cat$gAmpl[!is.na(cat$gamp)] <- cat$gamp[!is.na(cat$gamp)]

## boolean for rrab and rr
isrr <- !is.na(cat$subtype)
isab <- cat$subtype=="ab"
isab[is.na(isab)] <- FALSE


cols <- c(alpha("black",.01),alpha("red",1))

cex.lab <- 2
cex.axis <- 2
cex.legend <- 2

## use rr lyrae features estimated by sesar, tighter
png("../figs/sdss_color_period.png",height=700,width=800)
par(mar=c(5,5,1,1))
plot(cat$P,cat$ug,col=cols[isab+1],ylim=c(-.5,3.5),log="x",
     xlab="Period Estimate",ylab="ug color",cex.lab=cex.lab,pch=isab+1,axes=FALSE,
     xlim=c(0.2,5000),xaxs="i")
box()
axis(1,at=c(0.1,1,10,100,1000),labels=c(0.1,1,10,100,1000),cex.axis=cex.axis)
axis(2,cex.axis=cex.axis)
legend("topright",c("RR Lyrae","Not RR Lyrae"),col=2:1,pch=2:1,cex=cex.legend)
dev.off()


png("../figs/sdss_gamp_period.png",height=700,width=800)
par(mar=c(5,5,1,1))
plot(cat$P,cat$gAmpl,col=cols[isab+1],log="x",
     xlab="Period Estimate",ylab="g amplitude (peak to peak)",cex.lab=cex.lab,pch=isab+1,axes=FALSE,
     xlim=c(min(cat$P),1000),xaxs="i",ylim=c(0,2),yaxs='i')
box()
axis(1,at=c(0.1,1,10,100,1000),labels=c(0.1,1,10,100,1000),cex.axis=cex.axis)
axis(2,cex.axis=cex.axis)
legend("topright",c("RR Lyrae","Not RR Lyrae"),col=2:1,pch=2:1,cex=cex.legend)
dev.off()





