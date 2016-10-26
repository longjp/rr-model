rm(list=ls())
set.seed(1234)
library(tree)
library(tables)
library(MASS)
source("plot-partition.R")


generate_data <- function(n){
    nc1 <- rbinom(1,n,.5)
    nc2 <- n - nc1
    x <- matrix(rnorm(2*nc1,mean=0,sd=.1),nrow=nc1)
    shift <- runif(nc1)
    x <- cbind(1,x + cbind(shift,shift^2))
    nc3 <- rbinom(1,nc2,.6)
    nc2 <- nc2 - nc3
    y1 <- cbind(rnorm(nc3,mean=.2,sd=.1),rnorm(nc3,mean=.6,sd=.1))
    y2 <- cbind(rnorm(nc2,mean=.9,sd=.2),rnorm(nc2,mean=0,sd=.2))
    y <- cbind(2,rbind(y1,y2))
    dat <- rbind(x,y)
    colnames(dat) <- c("class","x1","x2")
    return(as.data.frame(dat))
}

dat <- generate_data(1000)
plot(dat[,2:3],col=dat[,1])


n <- 1000
dat <- generate_data(n)
test <- generate_data(100)


## plot training data
lims <- c(range(dat$x1),range(dat$x2))
pdf(paste0("../figs/training.pdf"))
par(mar=c(5,5,1,1))
plot(0,0,xlim=c(-.5,1.5),ylim=c(-.5,1.5),col=0,xlab="Feature x1",ylab="Feature x2",cex.lab=1.3)
points(dat$x1,dat$x2,col=as.character(dat$class),pch=as.numeric(dat$class))
legend("topleft",c("Class 1","Class 2"),col=1:2,pch=1:2,cex=1.5)
dev.off()




########### 1) Kernel density estimator classifier

## look at different levels of smoothing
hs <- c(0.04,.4,1.8)
for(ii in 1:length(hs)){
    pdf(paste0("../figs/kde2d",ii,".pdf"))
    par(mar=c(5,5,1,1))
    plot(0,0,xlim=c(-.5,1.5),ylim=c(-.5,1.5),col=0,xlab="Feature x1",ylab="Feature x2",cex.lab=1.3)
    points(dat$x1,dat$x2,col=as.character(dat$class),pch=as.numeric(dat$class))
    for(jj in 1:2){
        to_use <- dat$class==jj
        fit.kde <- kde2d(dat$x1[to_use],dat$x2[to_use],n=50,h=hs[ii],lims=c(-.5,1.5,-.5,1.5))
        contour(fit.kde,levels  =  seq(from=.2,to=.8,by=.3),add=TRUE,col=jj,cex=2)
    }
    legend("topleft",c("Class 1","Class 2"),col=1:2,pch=1:2,cex=1.5)
    dev.off()
}






### run cross validation on training


## computing density at new points is a pain, right now
## using strategy discussed here
## http://stackoverflow.com/questions/16201906/how-can-i-get-the-value-of-a-kernel-density-estimate-at-specific-points
compute_density <- function(x_old,y_old,x_new,y_new,h){
    dens <- kde2d(x_old,y_old,h=h,n=50)
    gr <- data.frame(with(dens, expand.grid(dens$x,dens$y)), as.vector(dens$z))
    names(gr) <- c("xgr", "ygr", "zgr")
    mod <- loess(zgr~xgr*ygr, data=gr)
    pred <- predict(mod, newdata=data.frame(xgr=x_new,ygr=y_new))
    pred[is.na(pred)] <- 0
    return(pred)
}

hs <- 0.0001*2^(0:20)
cv_group <- (1:nrow(dat))%%10 + 1
cv_error <- matrix(0,ncol=length(hs),nrow=10)

for(ii in 1:10){
    print(ii)
    for(jj in 1:length(hs)){
        x <- dat$x1[cv_group!=ii & dat$class==1]
        y <- dat$x2[cv_group!=ii & dat$class==1]
        preds1 <- compute_density(x,y,dat$x1[cv_group==ii],
                                  dat$x2[cv_group==ii],h=hs[jj])
        x <- dat$x1[cv_group!=ii & dat$class==2]
        y <- dat$x2[cv_group!=ii & dat$class==2]
        preds2 <- compute_density(x,y,dat$x1[cv_group==ii],
                                  dat$x2[cv_group==ii],h=hs[jj])
        cl.pred <- 1 + (preds2 > preds1)
        cv_error[ii,jj] <- mean(cl.pred != dat$class[cv_group==ii])
    }
}


## plot CV error as a function of bandwidth
pdf("../figs/kde_cv.pdf",width=10,height=6)
par(mar=c(5,5,1,1))
plot(.5,0,col=0,xlim=range(hs),ylim=c(0,max(cv_error)),
     xlab="Bandwidth",ylab="Missclassification Rate",log="x",
     cex.lab=1.3,xaxs='i')
matlines(hs,t(cv_error),col='black',lty=1)
points(hs,colMeans(cv_error),type='l',lwd=3,col='red')
abline(v=hs[which.min(colMeans(cv_error))],lwd=2)
bw_best <- hs[which.min(colMeans(cv_error))]
legend("topright",c("Cross Validation Run","Mean of Cross Validation"),col=1:2,lwd=2,cex=1.5)
dev.off()



## make predictions
pred1 <- compute_density(dat$x1[dat$class==1],dat$x2[dat$class==1],test$x1,test$x2,h=bw_best)
pred2 <- compute_density(dat$x1[dat$class==2],dat$x2[dat$class==2],test$x1,test$x2,h=bw_best)
pred <- (pred2 > pred1) + 1
table(pred,test$cl)

tblr <- tabular((Truth=as.factor(c("black","red")[test$class])) ~ (Predicted = as.factor(c("black","red")[pred])))
sink("../figs/confusion_kde.tex")
latex(tblr)
sink()








########### 2) CART classifier
dat$class <- c("black","red")[dat$class]
dat$class <- as.factor(dat$class)

fit.tree <- tree(class~.,data=dat)



pdf("../figs/sim_features.pdf",width=6,height=6)
par(mar=c(4.5,4.5,1,1))
plot(dat$x1,dat$x2,col=as.character(dat$class),xlab="Feature x",ylab="Feature y",cex.lab=1.5,
     pch=as.numeric(dat$class))
dev.off()

source('plot-partition.R')




## plot sequence of trees
n.best <- max(prune.tree(fit.tree)$size)
for(ii in 2:n.best){
    fit <- prune.tree(fit.tree,best=ii)
    pdf(paste("../figs/tree_",ii,".pdf",sep=""),width=12,height=6)
    par(mfcol=c(1,2),mar=c(4.5,4.5,1,1))
    plot(dat$x1,dat$x2,col=as.character(dat$class),xlab="Feature x1",ylab="Feature x2",cex.lab=1.5,
         pch=as.numeric(dat$class))
    partition.tree(fit,cex=1.5,add=TRUE,ordvars=c("x1","x2"),font=2)
    plot(fit,cex.lab=2,type="uniform")
    text(fit,digits=1,cex=1.5,font=2)
    dev.off()
}



##
library(rpart)
fit.tree <- rpart(class~.,data=dat,control=rpart.control(cp=.00001,minsplit=1))
pdf("../figs/cart_cv.pdf",width=10,height=8)
par(mar=c(5,5,5,1))
plotcp(fit.tree,cex.lab=1.3,cex.main=2)
dev.off()



