### plot outputs from ls periodogram applied to well and poorly sampled data
rm(list=ls())
library(multiband)
load("results.RData")
load("../../data/clean/sdss_sim_class.RData")


period_est_lomb_FULL <- period_est_lomb_FULL[,1]
period_est_lomb <- period_est_lomb[,1]

plot(period_est_lomb_FULL[cl=="rr"],period_est_lomb[cl=="rr"])
#### TODO: fix this code

## 1. determine N x p matrix for FULL
out_FULL <- array(0,dim=c(length(period_est_lomb_FULL),5,3),
             dimnames=list(names(period_est_lomb_FULL),
                           c("g","i","r","u","z"),
                           c("beta0","amp","rho")))
for(ii in 1:dim(out_FULL)[1]){
    tm <- tms_FULL[[ii]]
    p <- period_est_lomb_FULL[ii]
    test <- pgls(tm,periods=p,BCD_flag=FALSE)$best_fitLS
    out_FULL[ii,,] <- test
}


## 2. determine N x p matrix for downsampled
out <- array(0,dim=c(length(period_est_lomb_FULL),4,3),
             dimnames=list(names(period_est_lomb_FULL),
                           c("g","i","r","z"),
                           c("beta0","amp","rho")))
for(ii in 1:dim(out)[1]){
    tm <- tms[[ii]]
    p <- period_est_lomb[ii]
    l <- pgls(tm,periods=p,BCD_flag=FALSE)$best_fitLS
    for(jj in dimnames(out)[[2]]){
        if(jj %in% rownames(l))
            out[ii,jj,] <- l[jj,]
    }
}





## find features which separate RRL well
rands <- sample(1:dim(out)[1])
## xlim <- range(out_FULL[,"g","amp"])
## ylim <- range(out_FULL[,"i","beta0"] - out_FULL[,"g","beta0"])
mark <- 1:2
names(mark) <- c("not","rr")

## pdf("features_lomb_full.pdf",width=6,height=5)
## par(mar=c(5,5,1,1))
## plot(out_FULL[rands,"g","amp"],
##      out_FULL[rands,"i","beta0"] - out_FULL[rands,"g","beta0"],
##      col=mark[cl[rands]],pch=mark[cl[rands]],
##      xlim=xlim,ylim=ylim,
##      xlab="g amplitude",
##      ylab="i-g color",cex.lab=1.5)
## dev.off()

## pdf("features_lomb_downsampled.pdf",width=6,height=5)
## par(mar=c(5,5,1,1))
## plot(out[rands,"g","amp"],
##      out[rands,"i","beta0"] - out[rands,"g","beta0"],
##      col=mark[cl[rands]],pch=mark[cl[rands]],
##      xlim=xlim,ylim=ylim,
##      xlab="g amplitude",
##      ylab="i-g color",cex.lab=1.3)
## dev.off()




pdf("features_lomb_full.pdf",width=6,height=5)
par(mar=c(5,5,1,1))
plot(out_FULL[,"g","amp"],
     period_est_lomb_FULL,
     col=mark[cl],pch=mark[cl],
     xlab="g amplitude",
     ylab="period (days)",cex.lab=1.3,
     log='x',xlim=c(.01,5),ylim=c(.2,1))
legend("topright",c("RR Lyrae","Not RR Lyrae"),col=2:1,pch=2:1,cex=1.5)
dev.off()


pdf("features_lomb_downsampled.pdf",width=6,height=5)
par(mar=c(5,5,1,1))
plot(out[,"g","amp"],
     period_est_lomb,
     col=mark[cl],pch=mark[cl],
     cex.lab=1.3,
     xlab="g amplitude",
     ylab="period (days)",cex.lab=1.3,
     log='x',xlim=c(.01,5),ylim=c(.2,1))
dev.off()
