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
feat <- "E.B.V."
class_names <- 1:2
names(class_names) <- c("rr","not")
train <- sample(c(TRUE,FALSE),nrow(dat),replace=TRUE)
train_cl <- cl[train]
train_cl <- class_names[train_cl]
train_x <- dat[train,feat]
test_cl <- cl[!train]
test_cl <- class_names[test_cl]
test_x <- dat[!train,feat]

## estimate probs using standard technique
probs_est <- OneDKDEClassifier(train_cl,train_x,test_x)
## estimate probs using np mm
fb <- density(train_x[train_cl==1])
fb <- CreateApproxFun(fb)
FbD <- ecdf(train_x[train_cl==1])
params <- EstimateParams(test_x,FbD)
num <- (1-params$alpha_hat)*fb(test_x)
den <- num + params$alpha_hat*params$fs(test_x)
probs_mm <-  num / den

## compare results
## really! strange pattern
plot(probs_est,probs_mm)
## make roc curves for probs_true,probs_est,probs_mm
## fb = class 1 = positive
## fs = class 2 = negative
roc_est <- ComputeROCCurve(test_cl-1,probs_est)
roc_mm <- ComputeROCCurve(test_cl-1,probs_mm)

plot(0,0,col=0,type='l',
     xlab="False Positive Rate",ylab="True Positive Rate",
     xlim=c(0,1),ylim=c(0,1))
points(roc_est,type='l',lwd=2,col='blue')
points(roc_mm,type='l',lwd=2,col='orange')

## mixture model is classifying some very obviously negative (non RRL)
## as positive (ROC curve is fairly horizontal initiall), why



### build random forest on training, apply to test
dat_train <- dat[train,]
dat_test <- dat[!train,]
rf.fit <- randomForest(cl~.,data=dat_train)
probs_rf <- predict(rf.fit,dat_test,type='prob')[,2]

pairs(cbind(probs_rf,probs_est,probs_mm))
roc_rf <- ComputeROCCurve(test_cl-1,probs_rf)
plot(0,0,col=0,type='l',
     xlab="False Positive Rate",ylab="True Positive Rate",
     xlim=c(0,1),ylim=c(0,1))
points(roc_est,type='l',lwd=2,col='blue')
points(roc_mm,type='l',lwd=2,col='orange')
points(roc_rf,type='l',lwd=2,col='black')


## double check probs_mm calculation, appears messed up
## especially for phi
