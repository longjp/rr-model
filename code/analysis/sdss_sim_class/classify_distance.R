rm(list=ls())
set.seed(1234)

## load necessary libraries
library('parallel')
library('multiband')
library('randomForest')
load("../../fit_template/template.RData")
source("../../fit_template/fit_template.R")
source("../../common/funcs.R")
source("../funcs.R")


## data source
load("../../data/clean/sdss_sim_class.RData")
load("results.RData")
source("../params.R")

fig.dir <- "figs_classify"
unlink(fig.dir,recursive=TRUE)
dir.create(fig.dir)


period_est <- period_est[,1] ## just use best fit period

sigs <- unlist(lapply(tms,function(tm){lapply(tm,function(x){x$sigma})}))
hist(sigs[sigs<.1])

rss.n <- rep(0,N)
rss.n2 <- rep(0,N)
rss.n3 <- rep(0,N)
coeffs <- matrix(0,ncol=4,nrow=N)
cols <- matrix(0,ncol=2,nrow=N)
for(ii in 1:N){
    tm <- tms[[ii]]
    omega <- 1/period_est[ii]
    lc <- TMtoLC(tm)
    coes <- ComputeCoeffs(lc,omega,tem)
    coeffs[ii,] <- coes
    dev <- rep(0,length(tm))
    dev2 <- rep(0,length(tm))
    dev3 <- rep(0,length(tm))
    for(jj in 1:length(tm)){
        pred <- (coes[1] + tem$betas[jj] + coes[2]*tem$dust[jj]
            + coes[3]*tem$template_funcs[[jj]]((tm[[jj]][,1]*omega + coes[4]) %% 1))
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

features <- cbind(coeffs,period_est,rss.n,rss.n2,rss.n3,cols)
colnames(features) <- c("mu","E[B-V]","a","phi","period","rss.n","rss.n2","rss.n3","ig","ri")

features_sdss <- features
cl_sdss <- cl
save(features_sdss,cl_sdss,file="feats.RData")


## d1 <- density(features[cl=="rr",6],bw="SJ")
## d2 <- density(features[cl=="not",6],bw="SJ")
## xlim <- range(c(d1$x,d2$x))
## ylim <- range(c(d1$y,d2$y))
## plot(0,0,xlim=xlim,ylim=ylim,xlab="RSS",ylab="Density")
## points(d1$x,d1$y,col="black",type='l')
## points(d2$x,d2$y,col="red",type='l')
## legend("topright",c("rr","not rr"),lty=1,col=1:2)

cols <- c("#00000030",'red')
names(cols) <- c("not","rr")
pchs <- c(1,2)
names(pchs) <- c("not","rr")
pdf("RRfeatures.pdf")
pairs(features,col=cols[cl],pch=pchs[cl])
dev.off()


cols <- c(1,2)
rands <- sample(nrow(features))
names(cols) <- c("not","rr")

colnames(features)
pdf("features_template_downsampled.pdf",width=8,height=8)
par(mar=c(5,5,1,1))
plot(features[rands,5],features[rands,6],col=cols[cl[rands]],pch=pchs[cl[rands]],
     xlab="Period Estimate",ylab="RSS/n",cex.lab=1.3)
dev.off()


cl.plot <- as.numeric(as.factor(cl))
pdf(paste0(fig.dir,"/scatterplot_features.pdf"),height=6,width=7)
pairs(features,col=cl.plot,lab=c("Mean","E[B-V]","Amp.","Phase","Period","RSSn"),
      pch=cl.plot)
dev.off()

pdf(paste0(fig.dir,"/a_vs_dev.pdf"),height=6,width=7)
par(mar=c(5,5,1,1))
plot(features[,3],features[,6],col=cl.plot,pch=cl.plot,
     xlab="Amplitude",ylab="Model Residual",cex.lab=1.3)
dev.off()


pairs(cbind("amp"=log10(features[,"a"]),
            "period"=features[,"period"],
            "dev"=log10(features[,"rss.n"])),
      col=cl.plot,pch=cl.plot)

##pdf(paste0(fig.dir,"/a_vs_dev_log.pdf"),height=6,width=7)
par(mar=c(5,5,1,1))
plot(features[,3],features[,6],col=cl.plot,pch=cl.plot,
     xlab="Amplitude",ylab="RSS/n",cex.lab=1.3,
     log="xy",xlim=c(.01,5))
legend("topleft",c("RR Lyrae","Not RR Lyrae"),col=2:1,pch=2:1,cex=1.5)
dev.new()
par(mar=c(5,5,1,1))
plot(features[,3],features[,7],col=cl.plot,pch=cl.plot,
     xlab="Amplitude",ylab="RSS2/n",cex.lab=1.3,
     log="xy",xlim=c(.01,5))
legend("topleft",c("RR Lyrae","Not RR Lyrae"),col=2:1,pch=2:1,cex=1.5)
dev.new()
par(mar=c(5,5,1,1))
plot(features[,3],features[,8],col=cl.plot,pch=cl.plot,
     xlab="Amplitude",ylab="RSS3/n",cex.lab=1.3,
     log="xy",xlim=c(.01,5))
legend("topleft",c("RR Lyrae","Not RR Lyrae"),col=2:1,pch=2:1,cex=1.5)

##dev.off()


pdf(paste0(fig.dir,"/period_vs_EBV.pdf"),height=6,width=7)
par(mar=c(5,5,1,1))
plot(features[,5],features[,2],col=cl.plot,pch=cl.plot,
     xlab="Period Estimate",ylab="E[B-V] Estimate",cex.lab=1.3,
     ylim=c(-.5,2))
legend("topright",c("RR Lyrae","Not RR Lyrae"),col=2:1,pch=2:1,cex=1.5)
dev.off()

load("../../fit_template/template.RData")

ampg <- diff(range(tem$templates['g',]))/2

pdf("features_template_downsampled.pdf",height=5,width=6)
par(mar=c(5,5,1,1))
plot(features[,'a']*ampg,period_est,col=cl.plot,pch=cl.plot,
     xlab="g amplitude",ylab="period (days)",cex.lab=1.3,
     log='x',xlim=c(.01,5),ylim=c(.2,1))
legend("topleft",c("RR Lyrae","Not RR Lyrae"),col=2:1,pch=2:1,cex=1.5)
dev.off()



########## USE RANDOM FOREST TO FIND RR LYRAE, OUTPUT CLASSIFICATIONS
########## AND DISTANCE ESTIMATES (for input into coords.R)


dat <- na.roughfix(data.frame(cl,features[,-1]))
save(dat,file="features.RData")
rf.fit <- randomForest(cl~.,data=dat)
dat$predict <- predict(rf.fit)
dat$mu <- features[,"mu"]
to_use <- dat$predict=="rr"
rf_rr <- data.frame(ra=ra[to_use],dec=dec[to_use],d=10^(dat$mu[to_use]/5 + 1)/1000)
save(rf_rr,file="rf_rr.RData")



table(predict(rf.fit),cl)



########## PLOT DISTANCES FOR LIGHTCURVES
pdf("distance_comparison.pdf",width=7,height=7)
to_use <- distance <= 120
par(mar=c(5,5,1,1))
x <- distance[cl=="rr" & to_use]
y <- 10^((features[cl=="rr" & to_use,"mu"])/5 + 1)/1000
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
dev.off()
