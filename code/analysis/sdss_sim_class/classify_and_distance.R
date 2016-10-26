rm(list=ls())

## load necessary libraries
library('parallel')
library('multiband')
library('rpart')
library('randomForest')
library(rpart.plot)
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


## shouldn't tot.dev be mean.dev?

period_est <- period_est[,1] ## just use best fit period

tot.dev <- rep(0,N)
for(ii in 1:N){
    tm <- tms[[ii]]
    omega <- 1/period_est[ii]
    lc <- TMtoLC(tm)
    coeffs <- ComputeCoeffs(lc,omega,tem)
    dev <- rep(0,length(tm))
    for(jj in 1:length(tm)){
        pred <- (coeffs[1] + tem$betas[jj] + coeffs[2]*tem$dust[jj]
            + coeffs[3]*tem$template_funcs[[jj]]((tm[[jj]][,1]*omega + coeffs[4]) %% 1))
        dev[jj] <- sum(abs((pred - tm[[jj]][,2])))
    }
    tot.dev[ii] <- sum(dev)
}

coeffs <- matrix(0,ncol=4,nrow=N)
for(ii in 1:N){
    tm <- tms[[ii]]
    omega <- 1/period_est[ii]
    coeffs[ii,] <- ComputeCoeffs(TMtoLC(tm),omega,tem)
}


features <- cbind(coeffs,period_est,tot.dev)
colnames(features) <- c("mu","E[B-V]","a","phi","period","dev")


cols <- c("#00000030",'red')
names(cols) <- c("not","rr")
pchs <- c(1,2)
names(pchs) <- c("not","rr")
pdf("RRfeatures.pdf")
pairs(features,col=cols[cl],pch=pchs[cl])
dev.off()


dat <- data.frame(cl,features[,-1])
rpart.fit <- rpart(cl~.,data=dat,control=rpart.control(xval=10,minsplit=1,cp=.001))
rpart.fit$cptable
mincp <- rpart.fit$cptable[which.min(rpart.fit$cptable[,"xerror"]),"CP"]


pdf("real_data_cp.pdf",width=10,height=8)
par(mar=c(5,5,5,1))
plotcp(rpart.fit,cex.lab=1.5)
dev.off()


pdf("rpart_fulltree.pdf")
par(mar=c(0,1,0,1))
prp(rpart.fit,extra=2,compress=FALSE,varlen=0)
dev.off()

pdf("rpart_pruned.pdf")
rpart.fit.pruned <- prune(rpart.fit, cp= mincp)
prp(rpart.fit.pruned,extra=2,compress=FALSE,varlen=0)
dev.off()

pdf("rpart_pruned_more.pdf")
rpart.fit.pruned <- prune(rpart.fit, cp= 0.023)
prp(rpart.fit.pruned,extra=2,compress=FALSE,varlen=0)
dev.off()











out <- table(predict(rpart.fit.pruned,type="class"),cl)



sum(diag(out))/sum(out)

##pdf(paste0(fig.dir,"/cart_tree.pdf"),height=6,width=7)


par(xpd = TRUE)
plot(rpart.fit,uniform=TRUE)
text(rpart.fit, use.n = TRUE)
##dev.off()

rf.fit <- randomForest(cl~.,data=dat)

## predict everything with random forest
dat$predict <- predict(rf.fit)



cl.plot <- as.numeric(as.factor(cl))
pdf(paste0(fig.dir,"/scatterplot_features.pdf"),height=6,width=7)
pairs(features,col=cl.plot,lab=c("Mean","E[B-V]","Amp.","Phase","Period","Model Res."),
      pch=cl.plot)
dev.off()

pdf(paste0(fig.dir,"/a_vs_dev.pdf"),height=6,width=7)
par(mar=c(5,5,1,1))
plot(features[,3],features[,6],col=cl.plot,pch=cl.plot,
     xlab="Amplitude",ylab="Model Residual",cex.lab=1.3)
dev.off()





pairs(cbind("amp"=log10(features[,"a"]),
            "period"=features[,"period"],
            "dev"=log10(features[,"dev"])),
      col=cl.plot,pch=cl.plot)

pdf(paste0(fig.dir,"/a_vs_dev_log.pdf"),height=6,width=7)
par(mar=c(5,5,1,1))
plot(features[,3],features[,6],col=cl.plot,pch=cl.plot,
     xlab="Amplitude",ylab="Model Residual",cex.lab=1.3,
     log="xy",xlim=c(.01,5))
legend("topleft",c("RR Lyrae","Not RR Lyrae"),col=2:1,pch=2:1,cex=1.5)
dev.off()


pdf(paste0(fig.dir,"/period_vs_EBV.pdf"),height=6,width=7)
par(mar=c(5,5,1,1))
plot(features[,5],features[,2],col=cl.plot,pch=cl.plot,
     xlab="Period Estimate",ylab="E[B-V] Estimate",cex.lab=1.3)
legend("topright",c("RR Lyrae","Not RR Lyrae"),col=2:1,pch=2:1,cex=1.5)
dev.off()









######### now compute distances

fig.dir <- "figs_distance"
unlink(fig.dir,recursive=TRUE)
dir.create(fig.dir)


rrlyrae <- read.table("../../data/raw/apj326724t3_mrt.txt",skip=30)
names(rrlyrae)[1:5] <- c("ID","ra","dec","ar","d")
rrlyrae$ID <- paste0("LC_",rrlyrae$ID,".dat")


###### get coefficients AND classifications for ALL light curves here
rr_model <- data.frame(names(tms),dat)
names(rr_model)[1] <- "ID"
out <- merge(rrlyrae,rr_model,all=TRUE)
out <- out[!is.na(out$cl),] ## only use observations which we have classes for


## compute distance based on model (out$mu)

## pdf("distance_comparison.pdf")
## par(mar=c(5,5,1,1))

to_use <- out$predict=='rr'
plot(out$d[to_use],10^(out$mu[to_use]/5 + 1)/1000,col=(1*(out$cl[to_use]=='not')+1),
     xlab="Ground Truth Distance in kpc (Sesar 2010)",
     ylab="Estimate in kpc",
     cex.lab=1.3)
abline(a=0,b=1)

not_rr <- to_use & out$cl!='rr'
points(rep(min(out$d[to_use],na.rm=TRUE),sum(not_rr)),10^(out$mu[not_rr]/5 + 1)/1000,
       col='red',pch=19)
##dev.off()


#### make side by side comparison maps of Sesar

## why do we have rr, dec, d for some objects (5) that are not actually rr lyrae
## what happens to rrlyrae c and d? is my classifier wrong if it labels these as rrlyrae?`

out$d
plot(-out$d*sin(2*pi*out$ra/360),out$d*cos(2*pi*out$ra/360))


## go from l (longitude) ,b (latitude) (,d (distance)
## "Cartesian galactocentric coordinate system"

## def of stripe 82: RA ~ 20h to R.A. ∼ 4h (so 8 hours or 360/3 = 120 degrees)
