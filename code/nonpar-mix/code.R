rm(list=ls())
library(randomForest)
source("funcs.R")
set.seed(1234)

######### SIMULATION
## simulate two classes
n <- 50
alpha <- 0.25
class_means <- c(0,3)
train_cl <- sample(1:2,n,replace=TRUE,prob=c(1-alpha,alpha))
train_x <- rnorm(n,mean=class_means[train_cl],sd=1)
plot(density(train_x,bw="SJ"))
test_cl <- sample(1:2,n,replace=TRUE,prob=c(1-alpha,alpha))
test_x <- rnorm(n,mean=class_means[test_cl],sd=1)
plot(density(test_x,bw="SJ"))
## calculate probs using standard kde classifier
probs_est <- OneDKDEClassifier(train_cl,train_x,test_x)
## calculate true probs
num_true <- (1-alpha)*dnorm(test_x,mean=class_means[1],sd=1)
den_true <- num_true + alpha*dnorm(test_x,mean=class_means[2],sd=1)
probs_true <- num_true / den_true
## calculate probs using mixture model
fb <- density(train_x[train_cl==1])
fb <- CreateApproxFun(fb)
FbD <- ecdf(train_x[train_cl==1])
params <- EstimateParams(test_x,FbD)
num <- (1-params$alpha_hat)*fb(test_x)
den <- num + params$alpha_hat*params$fs(test_x)
probs_mm <-  num / den
probs <- cbind(probs_true,probs_est,probs_mm)
pairs(probs)

## make roc curves for probs_true,probs_est,probs_mm
## fb = class 1 = positive
## fs = class 2 = negative
roc_true <- ComputeROCCurve(test_cl-1,probs_true)
roc_est <- ComputeROCCurve(test_cl-1,probs_est)
roc_mm <- ComputeROCCurve(test_cl-1,probs_mm)
plot(0,0,col=0,type='l',
     xlab="False Positive Rate",ylab="True Positive Rate",
     xlim=c(0,1),ylim=c(0,1))
points(roc_true,type='l',lwd=2)
points(roc_est,type='l',lwd=2,col='blue')
points(roc_mm,type='l',lwd=2,col='orange')





######### "REAL" DATA
load("features.RData")
cl <- as.character(dat$cl)
##feat <- "phi"
feats <- names(dat)[2:length(dat)]
alphas <- rep(0,length(feats))
names(alphas) <- feats
## select training and test data
class_names <- 1:2
names(class_names) <- c("rr","not")
train <- sample(c(TRUE,FALSE),nrow(dat),replace=TRUE)
train_cl <- cl[train]
train_cl <- class_names[train_cl]
test_cl <- cl[!train]
test_cl <- class_names[test_cl]
## roc curves
roc_est <- array(0,dim=c(length(feats),length(test_cl),2),dimnames=list(feats,NULL,NULL))
roc_mm <- array(0,dim=c(length(feats),length(test_cl),2),dimnames=list(feats,NULL,NULL))
## probs
probs_est <- matrix(0,nrow=length(feats),ncol=length(test_cl))
probs_mm <- matrix(0,nrow=length(feats),ncol=length(test_cl))
## build random forest on training, apply to test
dat_train <- dat[train,]
dat_test <- dat[!train,]
rf.fit <- randomForest(cl~.,data=dat_train)
probs_rf <- predict(rf.fit,dat_test,type='prob')[,2]
roc_rf <- ComputeROCCurve(test_cl-1,probs_rf)
## run one class classifiers on every feature
for(ii in 1:length(feats)){
    feat <- feats[ii]
    train_x <- dat[train,feat]
    test_x <- dat[!train,feat]
    ## estimate probs using standard technique
    probs_est[ii,] <- OneDKDEClassifier(train_cl,train_x,test_x)
    ## estimate probs using np mm
    fb <- density(train_x[train_cl==1])
    fb <- CreateApproxFun(fb)
    FbD <- ecdf(train_x[train_cl==1])
    params <- EstimateParams(test_x,FbD)
    alphas[ii] <- params$alpha_hat
    num <- (1-params$alpha_hat)*fb(test_x)
    den <- num + params$alpha_hat*params$fs(test_x)
    probs_mm[ii,] <-  num / den
    roc_est[ii,,] <- ComputeROCCurve(test_cl-1,probs_est[ii,])
    roc_mm[ii,,] <- ComputeROCCurve(test_cl-1,probs_mm[ii,])
    pdf(paste0(feat,".pdf"))
    plot(0,0,col=0,type='l',
         xlab="False Positive Rate",ylab="True Positive Rate",
         xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i',
         main=feat)
    points(roc_est[ii,,],type='l',lwd=2,col='blue')
    points(roc_mm[ii,,],type='l',lwd=2,col='orange')
    points(roc_rf,type='l',lwd=2,col='black')
    legend("bottomright",c("Random Forest","KDE","Mixture Model"),
           col=c("black","blue","orange"),lty=1,lwd=2)
    dev.off()
}



pdf("roc_two_features.pdf",width=8,height=6)
par(mar=c(5,5,1,1))
plot(0,0,col=0,type='l',
     xlab="False Positive Rate",ylab="True Positive Rate",
     xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i',cex.lab=1.3)
points(roc_est["E.B.V.",,],type='l',lwd=2,lty=1,col='blue')
points(roc_mm["E.B.V.",,],type='l',lwd=2,lty=2,col='blue')
points(roc_est["ig",,],type='l',lwd=2,lty=1,col='orange')
points(roc_mm["ig",,],type='l',lwd=2,lty=2,col='orange')
points(roc_rf,type='l',lwd=2,col='black')
legend_names <- c("Random Forest (All Features)","KDE (E[B-V])","Mixture Model (E[B-V])","KDE (i-g)","Mixture Model (i-g)")
legend("bottomright",legend_names,
       col=c("black","blue","blue","orange","orange"),
       lty=c(1,1,2,1,2),lwd=2)
dev.off()


alphas <- alphas[order(alphas)]

pdf("dotchart.pdf",width=8,height=6)
par(mar=c(5,1,1,1))
dotchart(alphas,xlab="alpha hat",cex=1.3)
dev.off()

########## ATTEMPT TO BUILD 2-D KDE CLASSIFIER, FAILED
## library(MASS)


## dat2 <- dat_train
## dat2 <- dat2[,c("cl","E.B.V.","ig")]
## names(dat2) <- c("class","x1","x2")
## dat2$class <- as.numeric(dat2$class)
## head(dat2)


## dat3 <- dat_test
## dat3 <- dat3[,c("cl","E.B.V.","ig")]
## names(dat3) <- c("class","x1","x2")
## dat3$class <- as.numeric(dat3$class)
## head(dat3)




## ########### 1) Kernel density estimator classifier

## ## look at different levels of smoothing
## xlim <- c(-1,2)
## ylim <- c(-3.5,.5)
## hs <- c(0.04,.4,1.8)
## for(ii in 1:length(hs)){
##     pdf(paste0("kde2d",ii,".pdf"))
##     par(mar=c(5,5,1,1))
##     plot(0,0,col=0,
##          xlab="Feature x1",ylab="Feature x2",cex.lab=1.3,
##          xlim=xlim,ylim=ylim)
##     rands <- sample(1:nrow(dat2))
##     points(dat2$x1[rands],dat2$x2[rands],col=as.character(dat2$class[rands]),pch=as.numeric(dat2$class[rands]))
##     for(jj in 1:2){
##         to_use <- dat2$class==jj
##         fit.kde <- kde2d(dat2$x1[to_use],dat2$x2[to_use],
##                          n=50,h=hs[ii],lims=c(xlim,ylim))
##         contour(fit.kde,levels=seq(from=.2,to=.8,by=.3),
##                 add=TRUE,col=jj,cex=2)
##     }
##     legend("topleft",c("Class 1","Class 2"),col=1:2,pch=1:2,cex=1.5)
##     dev.off()
## }






## ### run cross validation on training


## ## computing density at new points is a pain, right now
## ## using strategy discussed here
## ## http://stackoverflow.com/questions/16201906/how-can-i-get-the-value-of-a-kernel-density-estimate-at-specific-points
## compute_density <- function(x_old,y_old,x_new,y_new,h){
##     dens <- kde2d(x_old,y_old,h=h,n=50)
##     gr <- data.frame(with(dens, expand.grid(dens$x,dens$y)), as.vector(dens$z))
##     names(gr) <- c("xgr", "ygr", "zgr")
##     mod <- loess(zgr~xgr*ygr, data=gr)
##     pred <- predict(mod, newdata=data.frame(xgr=x_new,ygr=y_new))
##     pred[is.na(pred)] <- 0
##     pred[pred<0] <- 0
##     return(pred)
## }

## hs <- 0.0001*2^(0:20)
## cv_group <- (1:nrow(dat2))%%10 + 1
## cv_error <- matrix(0,ncol=length(hs),nrow=10)

## for(ii in 1:10){
##     print(ii)
##     for(jj in 1:length(hs)){
##         x <- dat2$x1[cv_group!=ii & dat2$class==1]
##         y <- dat2$x2[cv_group!=ii & dat2$class==1]
##         preds1 <- compute_density(x,y,dat2$x1[cv_group==ii],
##                                   dat2$x2[cv_group==ii],h=hs[jj])
##         x <- dat2$x1[cv_group!=ii & dat2$class==2]
##         y <- dat2$x2[cv_group!=ii & dat2$class==2]
##         preds2 <- compute_density(x,y,dat2$x1[cv_group==ii],
##                                   dat2$x2[cv_group==ii],h=hs[jj])
##         cl.pred <- 1 + (preds2 > preds1)
##         cv_error[ii,jj] <- mean(cl.pred != dat2$class[cv_group==ii])
##     }
## }


## ## plot CV error as a function of bandwidth
## pdf("kde_cv.pdf",width=10,height=6)
## par(mar=c(5,5,1,1))
## plot(.5,0,col=0,xlim=range(hs),ylim=c(0,max(cv_error)),
##      xlab="Bandwidth",ylab="Missclassification Rate",log="x",
##      cex.lab=1.3,xaxs='i')
## matlines(hs,t(cv_error),col='black',lty=1)
## points(hs,colMeans(cv_error),type='l',lwd=3,col='red')
## abline(v=hs[which.min(colMeans(cv_error))],lwd=2)
## bw_best <- hs[which.min(colMeans(cv_error))]
## legend("topright",c("Cross Validation Run","Mean of Cross Validation"),col=1:2,lwd=2,cex=1.5)
## dev.off()



## ## make predictions
## alpha <- mean(dat2$class==1)
## pred1 <- compute_density(dat2$x1[dat2$class==1],
##                          dat2$x2[dat2$class==1],
##                          dat3$x1,dat3$x2,h=bw_best)
## pred2 <- compute_density(dat2$x1[dat2$class==2],
##                          dat2$x2[dat2$class==2],
##                          dat3$x1,dat3$x2,h=bw_best)
## num <- pred2*(1-alpha)
## den <- num + pred1*alpha
## probs_kde_2 <- num / den
## hist(probs_kde_2)



## plot(probs_est[1,],probs_kde_2)

## plot(probs_rf,probs_kde_2)


## roc_kde_2 <- ComputeROCCurve(test_cl-1,probs_kde_2)
## plot(roc_kde_2)















## ##pdf("roc_two_features.pdf",width=8,height=6)
## par(mar=c(5,5,1,1))
## plot(0,0,col=0,type='l',
##      xlab="False Positive Rate",ylab="True Positive Rate",
##      xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i',cex.lab=1.3)
## points(roc_est["E.B.V.",,],type='l',lwd=2,lty=1,col='blue')
## points(roc_mm["E.B.V.",,],type='l',lwd=2,lty=2,col='blue')
## points(roc_est["ig",,],type='l',lwd=2,lty=1,col='orange')
## points(roc_mm["ig",,],type='l',lwd=2,lty=2,col='orange')
## points(roc_rf,type='l',lwd=2,col='black')
## points(roc_kde_2,lwd=2,type='l',col="red")
## legend_names <- c("Random Forest (All Features)","KDE (E[B-V])","Mixture Model (E[B-V])","KDE (i-g)","Mixture Model (i-g)")
## legend("bottomright",legend_names,
##        col=c("black","blue","blue","orange","orange"),
##        lty=c(1,1,2,1,2),lwd=2)
## ##dev.off()
