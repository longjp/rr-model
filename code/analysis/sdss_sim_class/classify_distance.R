rm(list=ls())

## load necessary libraries
library('parallel')
library('multiband')
library('randomForest')
load("../../fit_template/template.RData")
source("../../fit_template/template.R")
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

rss.n <- rep(0,N)
coeffs <- matrix(0,ncol=4,nrow=N)
for(ii in 1:N){
    tm <- tms[[ii]]
    omega <- 1/period_est[ii]
    lc <- TMtoLC(tm)
    coes <- ComputeCoeffs(lc,omega,tem)
    coeffs[ii,] <- coes
    dev <- rep(0,length(tm))
    for(jj in 1:length(tm)){
        pred <- (coes[1] + tem$betas[jj] + coes[2]*tem$dust[jj]
            + coes[3]*tem$template_funcs[[jj]]((tm[[jj]][,1]*omega + coes[4]) %% 1))
        dev[jj] <- sum(abs((pred - tm[[jj]][,2])))
    }
    rss.n[ii] <- sum(dev) / sum(vapply(tm,nrow,c(0)))
}

features <- cbind(coeffs,period_est,rss.n)
colnames(features) <- c("mu","E[B-V]","a","phi","period","rss.n")


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
par(mar=c(5,5,1,1))

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

pdf(paste0(fig.dir,"/a_vs_dev_log.pdf"),height=6,width=7)
par(mar=c(5,5,1,1))
plot(features[,3],features[,6],col=cl.plot,pch=cl.plot,
     xlab="Amplitude",ylab="RSS/n",cex.lab=1.3,
     log="xy",xlim=c(.01,5))
legend("topleft",c("RR Lyrae","Not RR Lyrae"),col=2:1,pch=2:1,cex=1.5)
dev.off()


pdf(paste0(fig.dir,"/period_vs_EBV.pdf"),height=6,width=7)
par(mar=c(5,5,1,1))
plot(features[,5],features[,2],col=cl.plot,pch=cl.plot,
     xlab="Period Estimate",ylab="E[B-V] Estimate",cex.lab=1.3)
legend("topright",c("RR Lyrae","Not RR Lyrae"),col=2:1,pch=2:1,cex=1.5)
dev.off()



########## USE RANDOM FOREST TO FIND RR LYRAE, OUTPUT CLASSIFICATIONS
########## AND DISTANCE ESTIMATES (for input into coords.R)


dat <- data.frame(cl,features[,-1])
rf.fit <- randomForest(cl~.,data=dat)
dat$predict <- predict(rf.fit)
dat$mu <- features[,"mu"]
to_use <- dat$predict=="rr"
rf_rr <- data.frame(ra=ra[to_use],dec=dec[to_use],d=10^(dat$mu[to_use]/5 + 1)/1000)
save(rf_rr,file="rf_rr.RData")

table(predict(rf.fit),cl)



########## PLOT DISTANCES FOR LIGHTCURVES
dat <- data.frame(cl,features)

## fig.dir <- "figs_distance"
## unlink(fig.dir,recursive=TRUE)
## dir.create(fig.dir)

rrlyrae <- read.table("../../data/raw/apj326724t3_mrt.txt",skip=30)
names(rrlyrae)[1:5] <- c("ID","ra","dec","ar","d")
rrlyrae$ID <- paste0("LC_",rrlyrae$ID,".dat")


nrow(dat[dat$cl=="rr",])
###### get coefficients AND classifications for ALL light curves here
rr_model <- data.frame(names(tms),dat)
names(rr_model)[1] <- "ID"
out <- merge(rrlyrae,rr_model,all=TRUE)
out <- merge(rrlyrae,rr_model)


plot(out$d,10^(out$mu/5 + 1)/1000,
     xlab="Sesar 2010 Distance",
     ylab="Estimate from Sparsely Sampled",
     cex.lab=1.8)


sum(rrlyrae[,1] %in% rr_model[,1])

## why are there NAs here
lim <- range(c(out[out$cl=="rr","d"],10^(out[out$cl=="rr","mu"]/5 + 1)/1000))



pdf("distance_comparison.pdf",width=7,height=7)
par(mar=c(5,5,1,1))
plot(out[out$cl=="rr","d"],10^(out[out$cl=="rr","mu"]/5 + 1)/1000,
     xlab="Sesar 2010 Distance",
     ylab="Estimate from Sparsely Sampled",
     xlim=lim,ylim=lim,cex.lab=1.8)
dev.off()
