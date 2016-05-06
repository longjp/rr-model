rm(list=ls())
source('funcs.R')
library(scales)
set.seed(1234)

folder <- "../../../code/"

f <- paste0(folder,"apj326724t2_mrt.txt")
rrlyrae <- read.table(f,skip=42)
f <- list.files(paste0(folder,"rrlyrae"))

## get first rrlyrae
tms <- read.table(paste(folder,"rrlyrae/",f[1],sep=""))
p_est <- rrlyrae[rrlyrae[,1] == sub("LC_","",sub(".dat","",f[1])),3]
bands <- unique(tms[,2])
bands <- as.character(bands)
bands <- bands[order(bands)]

## plot unfolded
cex <- 1.3
xlim <- range(tms[,1])
ylim <- rev(range(tms[,3]))
pdf("../figures/unfolded.pdf",width=8,height=4)
par(mar=c(4.5,4.5,1,1))
plot(0,0,ylim=ylim,xlim=xlim,col=0,xlab="Time",ylab="Magnitude",cex.lab=cex,xaxs='i')
segments(tms[,1],tms[,3] - tms[,4],tms[,1],tms[,3] + tms[,4])
for(ii in 1:length(bands)){
    sub <- tms[,2] == bands[ii]
    points(tms[sub,1],tms[sub,3],col=ii,pch=ii,cex=cex)
}
legend("bottomleft",paste(as.character(bands)," band",sep=""),pch=1:length(bands),col=1:length(bands),cex=cex)
dev.off()

## plot folded light curve
xlim <- c(0,1)
ylim <- rev(range(tms[,3]))
K <- 3
pdf("../figures/folded.pdf",width=8,height=4)
par(mar=c(4.5,4.5,1,1))
plot(0,0,ylim=ylim,xlim=xlim,col=0,xlab=paste("Phase (period = ",round(p_est,2)," days )",sep=""),
     ylab="Magnitude",cex.lab=cex,xaxs='i')
segments((tms[,1] %% p_est)/p_est,tms[,3] - tms[,4],(tms[,1] %% p_est)/p_est,tms[,3] + tms[,4])
for(ii in 1:length(bands)){
    sub <- tms[,2] == bands[ii]
    points((tms[sub,1] %% p_est)/p_est,tms[sub,3],col=ii,pch=ii,cex=cex)
    ## draw model on plot
    lc <- cbind(tms[sub,1],tms[sub,3],1)
    X <- construct_design(2*pi/p_est,K,lc[,1])
    beta <- compute_params(2*pi/p_est,K,lc[,2],lc[,3]^2,X)
    N_t <- 100
    times <- seq(from=0,to=p_est,length.out=N_t)
    X <- construct_design(2*pi/p_est,K,times)
    points(times/p_est,X%*%beta,type='l',col=alpha(ii,.5))
}
dev.off()













########
######## panstarrs light curve
f <- list.files(paste0(folder,"PS1_sample_LCs"))

## get first rrlyrae
tms <- read.table(paste0(folder,"PS1_sample_LCs/",f[1],sep=""),header=TRUE)
period <- rrlyrae[rrlyrae[,1] == sub("LC_","",sub(".dat","",f[1])),3]
bands <- unique(tms[,2])
bands <- as.character(bands)
bands <- bands[order(bands)]

## plot unfolded
cex <- 1.3
xlim <- range(tms[,1])
ylim <- rev(range(tms[,3]))
pdf("../figures/unfolded_panstarrs.pdf",width=8,height=4)
par(mar=c(4.5,4.5,1,1))
plot(0,0,ylim=ylim,xlim=xlim,col=0,xlab="Time",ylab="Magnitude",cex.lab=cex,xaxs='i')
segments(tms[,1],tms[,3] - tms[,4],tms[,1],tms[,3] + tms[,4])
for(ii in 1:length(bands)){
    sub <- tms[,2] == bands[ii]
    points(tms[sub,1],tms[sub,3],col=ii,pch=ii,cex=cex)
}
legend("bottomleft",paste(as.character(bands)," band",sep=""),pch=1:length(bands),col=1:length(bands),cex=cex)
dev.off()
