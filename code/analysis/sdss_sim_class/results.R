rm(list=ls())

## load necessary libraries
library('parallel')
library('multiband')
load("../../fit_template/template.RData")
source("../../fit_template/template.R")
source("../../common/funcs.R")
source("../funcs.R")


## data source
load("../../data/clean/sdss_sim_class.RData")
load("results.RData")
source("../params.R")

unlink("figs",recursive=TRUE)
dir.create("figs")




## fraction of times true period in top 5
N <- sum(cl=="rr")
print(paste0("accuracies, top ",topN,":"))
print(paste("1%:",mean(within_x(period_est[1:N,],periods[1:N],0.01))))
print(paste("0.1%:",mean(within_x(period_est[1:N,],periods[1:N],0.001))))
print(paste("0.01%:",mean(within_x(period_est[1:N,],periods[1:N],0.0001))))
print("")

print(paste0("accuracies, top ",topN," for lomb:"))
print(paste("1%:",mean(within_x(period_est_lomb[1:N,],periods[1:N],0.01))))
print(paste("0.1%:",mean(within_x(period_est_lomb[1:N,],periods[1:N],0.001))))
print(paste("0.01%:",mean(within_x(period_est_lomb[1:N,],periods[1:N],0.0001))))
print("")

## fraction of times period is best
print("accuracies, top period:")
period_est <- period_est[,1] ## just use best fit period
print(paste("1%:",mean(abs((periods[1:N] - period_est[1:N])/periods[1:N]) < 0.01)))
print(paste("0.1%:",mean(abs((periods[1:N] - period_est[1:N])/periods[1:N]) < 0.001)))
print(paste("0.01%:",mean(abs((periods[1:N] - period_est[1:N])/periods[1:N]) < 0.0001)))

## fraction of times period is best
print("accuracies, top period, lomb:")
period_est_lomb <- period_est_lomb[,1]  ## just use best fit period
print(paste("1%:",mean(abs((periods[1:N] - period_est_lomb[1:N])/periods[1:N]) < 0.01)))
print(paste("0.1%:",mean(abs((periods[1:N] - period_est_lomb[1:N])/periods[1:N]) < 0.001)))
print(paste("0.01%:",mean(abs((periods[1:N] - period_est_lomb[1:N])/periods[1:N]) < 0.0001)))


lim <- range(c(periods[1:N],period_est_lomb[1:N],period_est[1:N]))
pdf("figs/period_comparison.pdf",height=8,width=18)
par(mfcol=c(1,2),mar=c(5,5,1,1))
plot(periods[1:N],period_est_lomb[1:N],xlim=lim,ylim=lim,
     xlab="True Period",ylab="Sine Period Estimate (Multiband Generalized Lomb Scargle)",
     cex.lab=1.5)
abline(a=0,b=1)
plot(periods[1:N],period_est[1:N],xlim=lim,ylim=lim,
     xlab="True Period",ylab="RR Lyrae Model Period Estimate",
     cex.lab=1.5)
abline(a=0,b=1)
dev.off()


## plot all bands with best fit parameters, store in figs
for(ii in 1:length(tms)){
    tm <- tms[[ii]]
    omega <- 1/period_est[ii]
    coeffs <- ComputeCoeffs(TMtoLC(tm),omega,tem)
    pdf(paste0("figs/",ii,".pdf"),height=12,width=8)
    par(mar=c(3,4,2,1),mfcol=c(5,1))
    for(jj in names(tem$dust)){
        pred <- (coeffs[1] + tem$betas[jj] + coeffs[2]*tem$dust[jj]
            + coeffs[3]*tem$template_funcs[[jj]]((tem$temp_time + coeffs[4]) %% 1))
        xlim <- range(tem$temp_time/omega)
        ylim <- range(range(pred),tm[[jj]]$mag)
        plot(tem$temp_time/omega,pred,type='l',xlab="Phase",ylab="Mag",ylim=rev(ylim),
         xlim=xlim,main=paste0(jj," band"))
        if(!is.null(tm[[jj]])){
            points((tm[[jj]]$time %% (1/omega)),tm[[jj]]$mag)
        }
    }
    dev.off()
}

## make all plots together
for(ii in 1:length(tms)){
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
