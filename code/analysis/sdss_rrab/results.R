rm(list=ls())

## load necessary libraries
## library('parallel')
## library('multiband')
library(RColorBrewer)
load("../../fit_template/template_sdss.RData")
source("../../fit_template/fit_template.R")
source("../../common/funcs.R")
source("../funcs.R")
source("../plot_funcs.R")

## data source
load("../../data/clean/sdss_rrab.RData")
source("../params.R")

unlink("figs",recursive=TRUE)
dir.create("figs")

coes <- matrix(0,nrow=length(tms),ncol=4)

## ## colors for plotting
## bandpch <- 1:6
## names(bandpch) <- c("u","g","r","i","z","Y")
## bandcol <- c("dodgerblue3","green","red",
##              "mediumorchid1","black","peachpuff4")
## names(bandcol) <- c("u","g","r","i","z","Y")


## plot all bands with best fit parameters, store in figs
rss.n <- vector("list",length(tms))
for(ii in 1:length(tms)){
    tm <- tms[[ii]]
    omega <- 1/periods[ii]
    coeffs <- ComputeCoeffs(TMtoLC(tm),omega,tem)
    coes[ii,] <- coeffs
    pdf(paste0("figs/",ii,".pdf"),height=12,width=8)
    par(mar=c(3,4,2,1),mfcol=c(5,1))
    dev <- vector("list",length(tm))
    for(jj in names(tem$dust)){
        pred <- (coeffs[1] + tem$betas[jj] + coeffs[2]*tem$dust[jj]
            + coeffs[3]*tem$template_funcs[[jj]]((tem$temp_time + coeffs[4]) %% 1))
        xlim <- range(tem$temp_time/omega)
        ylim <- range(range(pred),tm[[jj]]$mag)
        plot(tem$temp_time/omega,pred,type='l',xlab="Phase",ylab="Mag",ylim=rev(ylim),
         xlim=xlim,main=paste0(jj," band"))
        if(!is.null(tm[[jj]])){
            points((tm[[jj]]$time %% (1/omega)),tm[[jj]]$mag)
            segments((tm[[jj]]$time %% (1/omega)),tm[[jj]]$mag+tm[[jj]]$error,
            (tm[[jj]]$time %% (1/omega)),tm[[jj]]$mag-tm[[jj]]$error,
            col='grey')
        }
        pred <- (coeffs[1] + tem$betas[jj] + coeffs[2]*tem$dust[jj]
            + coeffs[3]*tem$template_funcs[[jj]]((tm[[jj]][,1]*omega + coeffs[4]) %% 1))
        dev[[jj]] <- abs(pred - tm[[jj]][,2])
    }
    rss.n[[ii]] <- unlist(dev)
    dev.off()
}


## make all plots together
for(ii in 1:length(tms)){
    print(ii)
    p_est <- periods[ii]
    omega <- 1 / p_est
    lc <- TMtoLC(tms[[ii]])
    coeffs <- ComputeCoeffs(lc,omega,tem)
    pdf(paste0("figs/",ii,"_one.pdf"),height=8,width=12)
    plotLC(lc,p_est,coeffs,tem)
    dev.off()
}






tab <- read.table("../../data/raw/apj326724t3_mrt.txt",skip=30)
tab <- tab[,c(1,4)]
names(tab) <- c("ID","dust")



dust <- data.frame(as.numeric(sub("LC_","",sub(".dat","",names(tms)))),coes[,2],coes[,1],rss.n)
names(dust) <- c("ID","my_dust","mu","rss")

nrow(tab)
nrow(dust)
out <- merge(tab,dust)
nrow(out)



## color points by distance




head(out)
lim <- range(c(out$dust,out$my_dust*tem$dust['r']))

pdf("dust_comparison.pdf",width=7,height=6)
par(mar=c(5,5,1,1))
plot(out$dust,out$my_dust*tem$dust['r'],
     xlim=lim,ylim=lim,
     xlab="Extinction in r-band from Schlegel",
     ylab="Extinction in r-band from RR Lyrae Model",cex.lab=1.3)
abline(a=0,b=1)
dev.off()


## make plots of best fitting and worst fitting light curves
ords <- order(rss.n,decreasing=FALSE)



## make all plots together
##for(ii in 1:length(tms)){
kk <- 0

kk <- kk + 1
ii <- ords[kk]
bands <- names(tem$dust)
bands <- sort(bands)
omega <- 1/periods[ii]
tm <- tms[[ii]]
lc <- TMtoLC(tm)
coeffs <- coes[ii,]
preds <- list()
for(jj in 1:length(bands)){
    preds[[jj]] <- (coeffs[1] + tem$betas[jj] + coeffs[2]*tem$dust[jj]
        + coeffs[3]*tem$template_funcs[[jj]]((tem$temp_time + coeffs[4]) %% 1))
}
ylim <- range(range(preds),range(lc$mag))
xlim <- range((lc$time %% (1/omega))*omega)
## pdf(paste0("figs/",ii,"_one.pdf"),height=8,width=12)
## par(mar=c(5,5,1,1))
plot(0,0,col=0,ylim=rev(ylim),xlim=xlim,xlab="Phase",ylab="Magnitude",cex.lab=1.5)
for(jj in 1:length(bands)){
    temp <- lc[lc$band==bands[jj],]
    if(nrow(temp) > 0.5){
        points((temp$time %% (1/omega)) * omega,temp$mag,col=jj,pch=jj,cex=1.5)
    }
    points(tem$temp_time,preds[[jj]],type='l',col=jj,lwd=3,lty=jj)
    segments((temp$time %% (1/omega)) * omega,temp$mag + 2*temp$error,
    (temp$time %% (1/omega)) * omega,temp$mag - 2*temp$error,col='grey')
}
legend("bottomleft",bands,col=1:length(bands),lty=1:length(bands),lwd=3,cex=1.5)

##    dev.off()
##}


## DO a weighted mean of residuals

