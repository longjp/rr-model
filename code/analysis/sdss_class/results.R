rm(list=ls())

## load necessary libraries
library('parallel')
library('multiband')
library('rpart')
library('randomForest')
load("../../fit_template/template_sdss.RData")
source("../../fit_template/fit_template.R")
source("../../common/funcs.R")
source("../funcs.R")


## data source
load("../../data/clean/sdss_class.RData")
load("results.RData")
source("../params.R")

fig.dir <- "figs_classify"
unlink(fig.dir,recursive=TRUE)
dir.create(fig.dir)


period_est <- period_est[,1] ## just use best fit period

tot.dev <- rep(0,N)
for(ii in 1:N){
    print(ii)
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
colnames(features) <- c("alpha","d","a","phi","period","dev")


cols <- 1:2
names(cols) <- c("rr","not")
pairs(features,col=cols[cl])

dat <- data.frame(cl,features)
rpart.fit <- rpart(cl~.,data=dat)
pdf(paste0(fig.dir,"/cart_tree.pdf"),height=6,width=7)
par(xpd = TRUE)
plot(rpart.fit,uniform=TRUE)
text(rpart.fit, use.n = TRUE)
dev.off()

rf.fit <- randomForest(cl~.,data=dat)


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


