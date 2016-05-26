rm(list=ls())

## load necessary libraries
library('parallel')
library('multiband')
load("../fit_template/template.RData")
source("../fit_template/template.R")
source("../common/funcs.R")
source("funcs.R")

## data source
load("../data/clean/sdss_sim.RData")
load("sdss_sim_period_est.RData")

unlink("figs",recursive=TRUE)
dir.create("figs")

## fraction of times true period in top 5
N <- nrow(period_est)
print("accuracies, top 5:")
print(paste("1%:",mean(within_x(period_est,periods[1:N],0.01))))
print(paste("0.1%:",mean(within_x(period_est,periods[1:N],0.001))))
print(paste("0.01%:",mean(within_x(period_est,periods[1:N],0.0001))))
print("")


N <- nrow(period_est)
print("accuracies, top 5 for lomb:")
print(paste("1%:",mean(within_x(period_est_lomb,periods[1:N],0.01))))
print(paste("0.1%:",mean(within_x(period_est_lomb,periods[1:N],0.001))))
print(paste("0.01%:",mean(within_x(period_est_lomb,periods[1:N],0.0001))))
print("")


## ## how often is algorithm (as opposed to model) is failing
## ## ie true period has lower rss than rss of top 5
## ## but true period is not in top 5
## rss_true <- rep(0,N)
## rss_est <- rep(0,N)
## phis <- (1:100)/100
## for(ii in 1:length(rss_true)){
##     tm <- tms[[ii]]
##     rss_est[ii] <- min(ComputeRSSPhase(tm,1/param$period[ii],tem,phis))
##     rss_true[ii] <- min(vapply(1/period_est[ii,],function(x){ComputeRSSPhase(tm,x,tem,phis)},rep(0,length(phis))))
## }
## is_correct <- within_x(period_est,periods[1:N],0.0001)
## print(paste0("fraction time wrong and correct period has lower rss: ",
##              round(mean((rss_true < rss_est) & !is_correct),4)))
## print("")

## fraction of times period is best
print("accuracies, top period:")
period_est <- period_est[,1] ## just use best fit period
print(paste("1%:",mean(abs((periods[1:N] - period_est)/periods[1:N]) < 0.01)))
print(paste("0.1%:",mean(abs((periods[1:N] - period_est)/periods[1:N]) < 0.001)))
print(paste("0.01%:",mean(abs((periods[1:N] - period_est)/periods[1:N]) < 0.0001)))

## fraction of times period is best
print("accuracies, top period:")
period_est_lomb <- period_est_lomb[,1]  ## just use best fit period
print(paste("1%:",mean(abs((periods[1:N] - period_est_lomb)/periods[1:N]) < 0.01)))
print(paste("0.1%:",mean(abs((periods[1:N] - period_est_lomb)/periods[1:N]) < 0.001)))
print(paste("0.01%:",mean(abs((periods[1:N] - period_est_lomb)/periods[1:N]) < 0.0001)))

## plot all bands with best fit parameters, store in figs
for(ii in 1:N){
    tm <- tms[[ii]]
    omega <- 1/period_est[ii]
    coeffs <- ComputeCoeffs(TMtoLC(tm),omega,tem)
    pdf(paste0("figs/",ii,".pdf"),height=12,width=8)
    par(mar=c(3,4,2,1),mfcol=c(5,1))
    for(jj in 1:length(tem$dust)){
        pred <- (coeffs[1] + tem$betas[jj] + coeffs[2]*tem$dust[jj]
            + coeffs[3]*tem$template_funcs[[jj]]((tem$temp_time + coeffs[4]) %% 1))
        plot(tem$temp_time/omega,pred,type='l',xlab="Phase",ylab="Mag",ylim=rev(range(pred)),main=paste0(names(tem$betas)[jj]," band"))
        points((tm[[jj]]$time %% (1/omega)),tm[[jj]]$mag)
    }
    dev.off()
}

## make all plots together
for(ii in 1:N){
    bands <- names(tem$dust)
    bands <- sort(bands)
    omega <- 1/period_est[ii]
    tm <- tms[[ii]]
    lc <- TMtoLC(tm)
    coeffs <- ComputeCoeffs(lc,omega,tem)
    preds <- list()
    for(jj in 1:length(bands)){
        preds[[jj]] <- (coeffs[1] + tem$betas[jj] + coeffs[2]*tem$dust[jj]
            + coeffs[3]*tem$template_funcs[[jj]]((tem$temp_time + coeffs[4]) %% 1))
    }
    ylim <- range(range(preds),range(lc$mag))
    xlim <- range((lc$time %% (1/omega))*omega)
    pdf(paste0("figs/",ii,"_one.pdf"),height=8,width=12)
    par(mar=c(5,5,1,1))
    plot(0,0,col=0,ylim=rev(ylim),xlim=xlim,xlab="Phase",ylab="Magnitude",cex.lab=1.5)
    for(jj in 1:length(bands)){
        temp <- lc[lc$band==bands[jj],]
        if(nrow(temp) > 0.5){
            points((temp$time %% (1/omega)) * omega,temp$mag,col=jj,pch=jj,cex=1.5)
        }
        points(tem$temp_time,preds[[jj]],type='l',col=jj,lwd=3,lty=jj)
    }
    legend("bottomleft",bands,col=1:length(bands),lty=1:length(bands),lwd=3,cex=1.5)
    dev.off()
}
