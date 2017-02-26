## M-x font-lock-mode

rm(list=ls())
library(scales)

fold <- "../../../code/data/raw/AllLCs/"
rrlyrae <- read.table("apj326724t2_mrt.txt",skip=42)
cat <- read.table("stripe82candidateVar_v1.1.dat",header=TRUE)
cat[,1] <- paste0("LC_",cat[,1],".dat")
fs <- list.files(fold)
fs <- fs[fs!="LC_reorganize.tcl.dat"] ## get rid of non lc file

cat <- cat[cat[,1] %in% fs,]
fs <- fs[fs %in% cat[,1]]
cat <- cat[order(cat[,1]),]
fs <- fs[order(fs)]

rrlyrae <- rrlyrae[paste0("LC_",rrlyrae[,1],".dat") %in% fs,]


ComputeColor <- function(lc,b1,b2){
    return(median(lc[lc[,2]==b1,3]) - median(lc[lc[,2]==b2,3]))
}

lcs <- vector("list",length(fs))
for(ii in 1:length(fs)){
    lcs[[ii]] <- read.table(paste0(fold,fs[ii]))
}

## compute colors
b1 <- "u"
b2 <- "g"
ug <- vapply(lcs,ComputeColor,c(0),b1=b1,b2=b2)
b1 <- "g"
b2 <- "r"
gr <- vapply(lcs,ComputeColor,c(0),b1=b1,b2=b2)
b1 <- "r"
b2 <- "i"
ri <- vapply(lcs,ComputeColor,c(0),b1=b1,b2=b2)
b1 <- "i"
b2 <- "z"
iz <- vapply(lcs,ComputeColor,c(0),b1=b1,b2=b2)



## use rr lyrae periods estimated by sesar, tighter
isrr <- cat$ID %in% paste0("LC_",rrlyrae$V1,".dat")
cat$rrlyrae <- 1*(isrr)
cat$P[isrr] <- rrlyrae[,3]


cols <- c(alpha("black",.01),alpha("red",1))
pdf("../figs/sdss_color_period.pdf",height=6,width=7)
par(mar=c(5,5,1,1))
plot(cat$P,ug,col=cols[isrr+1],ylim=c(-.5,3.5),log="x",
     xlab="Period Estimate",ylab="U-G (Color)",cex.lab=1.5,pch=isrr+1,axes=FALSE,
     xlim=c(0.1,5000),xaxs="i")
box()
axis(1,at=c(0.1,1,10,100,1000),labels=c(0.1,1,10,100,1000))
axis(2)
legend("topright",c("RR Lyrae","Not RR Lyrae"),col=2:1,pch=2:1,cex=1.5)
dev.off()
