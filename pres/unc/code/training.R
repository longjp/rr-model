rm(list=ls())
set.seed(1234)
library(rpart)
library(rpart.plot)
library(xtable)
library(tables)
library(tree)
n <- 100

## well-train
build_feats <- function(){
    rr.p <- rnorm(n,mean=0.5,sd=.15)
    rr.a <- runif(n,min=.2,max=.8)
    rr <- cbind(rr.p,rr.a)
    mira.p <- 100 + rexp(n,1/50)
    mira.a <- 1 + runif(n)
    mira <- cbind(mira.p,mira.a)
    feats <- rbind(rr,mira)
    feats <- data.frame(c(rep("rr",n),rep("mira",n)),feats)
    names(feats) <- c("class","period","amplitude")
    return(feats)
}

build_noisy <- function(){
    rr.p <- rexp(n,1/40)
    rr.a <- runif(n,min=0,max=1.2)
    rr <- cbind(rr.p,rr.a)
    mira.p <- rexp(n,1/40)
    mira.a <- 1 + runif(n,min=-.2,max=1)
    mira <- cbind(mira.p,mira.a)
    feats <- rbind(rr,mira)
    feats <- data.frame(c(rep("rr",n),rep("mira",n)),feats)
    names(feats) <- c("class","period","amplitude")
    return(feats)
}    


train1 <- build_feats()
train2 <- build_noisy()
unlabeled <- build_noisy()


## build tree on iid and clean training data
fit.rpart <- rpart(class~.,data=train1)
pdf("figs/train1_tree.pdf")
par(mar=c(0,0,0,0))
prp(fit.rpart,extra=2,compress=FALSE,varlen=0)
dev.off()

fit.rpart2 <- prune(rpart(class~.,data=train2),cp=.1)
pdf("figs/train2_tree.pdf")
prp(fit.rpart2,extra=2,compress=FALSE,varlen=0)
dev.off()


tblr <- tabular((Truth=unlabeled$class) ~ (Predicted = predict(fit.rpart,unlabeled,type="class")))
sink("figs/confusion_train1.tex")
latex(tblr)
sink()

tblr <- tabular((Truth=unlabeled$class) ~ (Predicted = predict(fit.rpart2,unlabeled,type="class")))
sink("figs/confusion_train2.tex")
latex(tblr)
sink()

## train 1 with split
xlim <- log(range(c(train1$period,train2$period,unlabeled$period)),base=10)
ylim <- range(c(train1$amplitude,train2$amplitude,unlabeled$amplitude))
pdf("figs/train1_features_vline.pdf",width=5,height=4)
par(mar=c(4,4,1,1))
plot(log(train1$period,base=10),train1$amplitude,col=c("blue","orange")[train1$class],
     pch=(1:2)[train1$class],xlab="log(period)",ylab="amplitude",cex.lab=1.2,xlim=xlim,ylim=ylim)
abline(v=log(51.06939,base=10),col="red",lwd=2)
legend("topleft",as.character(sort(unique(train1$class))),col=c("blue","orange"),pch=1:2,cex=1.2)
dev.off()

## unlabeled with split from train 1
xlim <- log(range(c(train1$period,train2$period,unlabeled$period)),base=10)
ylim <- range(c(train1$amplitude,train2$amplitude,unlabeled$amplitude))
pdf("figs/unlabeled_features_vline1.pdf",width=5,height=4)
par(mar=c(4,4,1,1))
plot(log(unlabeled$period,base=10),unlabeled$amplitude,col=c("blue","orange")[unlabeled$class],
     pch=(1:2)[unlabeled$class],xlab="log(period)",ylab="amplitude",cex.lab=1.2,xlim=xlim,ylim=ylim)
abline(v=log(51.06939,base=10),col="red",lwd=2)
legend("topleft",as.character(sort(unique(unlabeled$class))),col=c("blue","orange"),pch=1:2,cex=1.2)
dev.off()

## train 2 split
xlim <- log(range(c(train1$period,train2$period,unlabeled$period)),base=10)
ylim <- range(c(train1$amplitude,train2$amplitude,unlabeled$amplitude))
pdf("figs/train2_features_vline.pdf",width=5,height=4)
par(mar=c(4,4,1,1))
plot(log(train2$period,base=10),train2$amplitude,col=c("blue","orange")[train2$class],
     pch=(1:2)[train2$class],xlab="log(period)",ylab="amplitude",cex.lab=1.2,xlim=xlim,ylim=ylim)
abline(h=0.8333952,col="red",lwd=2)
legend("topleft",as.character(sort(unique(train2$class))),col=c("blue","orange"),pch=1:2,cex=1.2)
dev.off()

## unlabeled with train 2 split
xlim <- log(range(c(train1$period,train2$period,unlabeled$period)),base=10)
ylim <- range(c(train1$amplitude,train2$amplitude,unlabeled$amplitude))
pdf("figs/unlabeled_features_vline2.pdf",width=5,height=4)
par(mar=c(4,4,1,1))
plot(log(unlabeled$period,base=10),unlabeled$amplitude,col=c("blue","orange")[unlabeled$class],
     pch=(1:2)[unlabeled$class],xlab="log(period)",ylab="amplitude",cex.lab=1.2,xlim=xlim,ylim=ylim)
abline(h=0.8333952,col="red",lwd=2)
legend("topleft",as.character(sort(unique(unlabeled$class))),col=c("blue","orange"),pch=1:2,cex=1.2)
dev.off()


### make scatterplots of features
xlim <- log(range(c(train1$period,train2$period,unlabeled$period)),base=10)
ylim <- range(c(train1$amplitude,train2$amplitude,unlabeled$amplitude))
pdf("figs/train1_features.pdf",width=5,height=4)
par(mar=c(4,4,1,1))
plot(log(train1$period,base=10),train1$amplitude,col=c("blue","orange")[train1$class],
     pch=(1:2)[train1$class],xlab="log(period)",ylab="amplitude",cex.lab=1.2,xlim=xlim,ylim=ylim)
legend("topleft",as.character(sort(unique(train1$class))),col=c("blue","orange"),pch=1:2,cex=1.2)
dev.off()

pdf("figs/train2_features.pdf",width=5,height=4)
par(mar=c(4,4,1,1))
plot(log(train2$period,base=10),train2$amplitude,col=c("blue","orange")[train2$class],
     pch=(1:2)[train2$class],xlab="log(period)",ylab="amplitude",cex.lab=1.2,xlim=xlim,ylim=ylim)
legend("topleft",as.character(sort(unique(train2$class))),col=c("blue","orange"),pch=1:2,cex=1.2)
dev.off()

pdf("figs/unlabeled_features.pdf",width=5,height=4)
par(mar=c(4,4,1,1))
plot(log(unlabeled$period,base=10),unlabeled$amplitude,col=c("blue","orange")[unlabeled$class],
     pch=(1:2)[unlabeled$class],xlab="log(period)",ylab="amplitude",cex.lab=1.2,xlim=xlim,ylim=ylim)
legend("topleft",as.character(sort(unique(unlabeled$class))),col=c("blue","orange"),pch=1:2,cex=1.2)
dev.off()
