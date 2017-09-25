rm(list=ls())

## load necessary libraries
library('parallel')
library('multiband')
load("../../fit_template/template_sdss.RData")
source("../../fit_template/fit_template.R")
source("../../common/funcs.R")
source("../funcs.R")
source("../../common/plot_funcs.R")

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
pdf("period_est_sine.pdf",height=6,width=6)
par(mar=c(5,5,1,1))
plot(periods[1:N],period_est_lomb[1:N],xlim=lim,ylim=lim,
     xlab="True Period",ylab="Sine Period Estimate",
     cex.lab=1.5)
abline(a=0,b=1)
dev.off()

pdf("period_est_template.pdf",height=6,width=6)
par(mar=c(5,5,1,1))
plot(periods[1:N],period_est[1:N],xlim=lim,ylim=lim,
     xlab="True Period",ylab="RR Lyrae Model Period Estimate",
     cex.lab=1.5)
abline(a=0,b=1)
dev.off()



lim <- range(c(periods[1:N],period_est_lomb[1:N],period_est[1:N]))
pdf("period_est_comparison.pdf",height=6,width=12)
par(mar=c(5,5,1,1),mfcol=c(1,2))
plot(periods[1:N],period_est_lomb[1:N],xlim=lim,ylim=lim,
     xlab="True Period",ylab="Sine Period Estimate",
     cex.lab=2)
abline(a=0,b=1)
par(mar=c(5,5,1,1))
plot(periods[1:N],period_est[1:N],xlim=lim,ylim=lim,
     xlab="True Period",ylab="RRL Model Period Estimate",
     cex.lab=2)
abline(a=0,b=1)
dev.off()


## make all plots together
for(ii in 1:length(tms)){
    tm <- tms[[ii]]
    lc <- TMtoLC(tm)
    p_est <- period_est[ii]
    omega <- 1 / p_est
    coeffs <- ComputeCoeffs(lc,omega,tem)
    pdf(paste0("figs/",ii,"_one.pdf"),height=8,width=12)
    plotLC(lc,p_est,coeffs,tem,main=NULL,tem_only=TRUE)
    dev.off()
}

