## determine period estimation accuracy and
## distance estimation accuracy for full light curves


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



## fraction of times period is best
print("accuracies, top period:")
est <- period_est_FULL[,1] ## just use best fit period
print(paste("1%:",mean(abs((periods[1:N] - est[1:N])/periods[1:N]) < 0.01)))
print(paste("0.1%:",mean(abs((periods[1:N] - est[1:N])/periods[1:N]) < 0.001)))
print(paste("0.01%:",mean(abs((periods[1:N] - est[1:N])/periods[1:N]) < 0.0001)))


plot(periods[1:N],est[1:N],xlim=c(.2,1),ylim=c(.2,1))



###### estimate coefficients
N <- sum(cl=="rr")
rss.n <- rep(0,N)
rss.n2 <- rep(0,N)
rss.n3 <- rep(0,N)
coeffs <- matrix(0,ncol=4,nrow=N)
cols <- matrix(0,ncol=2,nrow=N)
for(ii in 1:N){
    tm <- tms[[ii]]
    omega <- 1/est[ii]
    lc <- TMtoLC(tm)
    coes <- ComputeCoeffs(lc,omega,tem)
    coeffs[ii,] <- coes
    dev <- rep(0,length(tm))
    dev2 <- rep(0,length(tm))
    dev3 <- rep(0,length(tm))
    for(jj in 1:length(tm)){
        pred <- PredictSingleBand(tm[[jj]][,1],names(tm)[jj],omega,coes,tem)
        ## pred <- (coes[1] + tem$betas[jj] + coes[2]*tem$dust[jj]
        ##     + coes[3]*tem$template_funcs[[jj]]((tm[[jj]][,1]*omega + coes[4]) %% 1))
        dev[jj] <- sum(abs((pred - tm[[jj]][,2])))
        dev2[jj] <- sum(abs((pred - tm[[jj]][,2])/tm[[jj]][,3]))
        dev3[jj] <- sum(abs((pred - tm[[jj]][,2])/(tm[[jj]][,3] + .2)))
    }
    rss.n[ii] <- sum(dev) / sum(vapply(tm,nrow,c(0)))
    rss.n2[ii] <- sum(dev2) / sum(vapply(tm,nrow,c(0)))
    rss.n3[ii] <- sum(dev3) / sum(vapply(tm,nrow,c(0)))
    mean_mags <- vapply(tm,function(x){mean(x[,2])},c(0))
    cols[ii,] <- c(mean_mags['i'] - mean_mags['g'],mean_mags['r'] - mean_mags['i'])
}
features <- cbind(coeffs,periods[1:N],rss.n,rss.n2,rss.n3,cols)
colnames(features) <- c("mu","E[B-V]","a","phi","period","rss.n","rss.n2","rss.n3","ig","ri")







to_use <- distance[1:N] <= 120
par(mar=c(5,5,1,1))
x <- distance[1:N][to_use]
y <- 10^((features[to_use,"mu"])/5 + 1)/1000
a <- features[to_use,"a"]
lim <- c(0,120)
plot(x,y,
     xlab="Sesar 2010 Distance (kpc)",
     ylab="Estimate from Sparsely Sampled (kpc)",
     cex.lab=1.8,xlim=lim,ylim=lim,
     xaxs='i',yaxs='i')
abline(a=0,b=1,lwd=2)
error_budget <- 0.1
lines(0:120,(1-error_budget)*(0:120),lty=2,lwd=2)
lines(0:120,(1+error_budget)*(0:120),lty=2,lwd=2)
legend("topleft",c("Identity",paste0(100*error_budget, "% Scatter")),lty=1:2,lwd=2,cex=1.5)
##dev.off()





## make all plots together
N <- sum(cl=="rr")
for(ii in 1:N){
    tm <- tms_FULL[[ii]]
    lc <- TMtoLC(tm)
    p_est <- periods[ii]
    omega <- 1 / p_est
    coeffs <- ComputeCoeffs(lc,omega,tem)
    pdf(paste0("figs/",ii,"_one.pdf"),height=8,width=12)
    plotLC(lc,p_est,coeffs,tem,main=NULL,tem_only=TRUE)
    dev.off()
}




head(features)
