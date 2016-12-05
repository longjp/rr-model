rm(list=ls())
source("funcs.R")
##load("features.RData")
set.seed(1234)

## TODO:
## compare oneD classifier to non-parametric mm on simulated
## compare 1 dimensional kde classifiers to non-parametric mm method
## compare 1 dimensional classifiers to random forest







## simulate two classes
n <- 1000
alpha <- 0.5
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



xo <- sort(x)

plot(xo,(1:n)/n)



cdf <- ((1:n)/n - (1-alpha)*pnorm(xo)) / alpha
out <- pava(cdf)
out[out<0] <- 0
out[out>1] <- 1

plot(xo,cdf,type='l',lwd=2)
points(xo,out,type='l',col='blue',lwd=2)
points(xo,pnorm(xo,mean=3),col="orange",lwd=2,type='l')


FbD <- pnorm
params <- EstimateParams(x,FbD)
probs <- (params$alpha_hat*params$fs(x)) / params$f(x)
hist(probs)


probs <- (params$alpha_hat*params$fs(x)) / (params$alpha_hat*params$fs(x) + (1-params$alpha_hat)*dnorm(x))


ComputeROCCurve <- function(cl,probs){
    cl_pred <- cl[order(probs,decreasing=TRUE)]
    return(cbind(cumsum(-cl_pred + 1)/sum(cl==0),cumsum(cl_pred)/sum(cl==1)))
}

plot(ComputeROCCurve(cl-1,probs),type='l')


probs_true <- alpha*dnorm(x,mean=3) / (alpha*dnorm(x,mean=3) + (1-alpha)*dnorm(x,mean=0))
plot(probs_true,probs)
abline(a=0,b=1)

par(mfcol=c(2,1))
plot(x,probs_true)
abline(v=1.5)
plot(x,probs)
abline(v=1.5)
abline(v=x[cl==1],col='grey')




roc_est <- ComputeROCCurve(cl-1,probs)
roc_true <- ComputeROCCurve(cl-1,probs_true)
plot(roc_est,type='l')
points(roc_true,col='red',type='l')

out <- EstMixMdl(x,FbD)
d <- out$distance
d2 <- ComputeSecondDeriv(d)
d2 <- d2 / max(d2)
d <- d / max(d)
plot(out$alpha_grid,d,type='l',xaxs='i')
points(out$alpha_grid,d2,lty=2,type='l')
abline(v=alpha)
alpha_hat <- out$alpha_grid[which.max(d2)]
alpha_hat
alpha


cdf_est <- out$fs[which.max(d2),]
plot(xo,cdf_est)
x2 <- xo[diff(c(0,cdf_est)) > 0]
d2 <- density(x2)
d2 <- approxfun(d2$x,d2$y)
plot((0:200)/33,d2((0:200)/33))
d <- density(xo)
d1 <- approxfun(d$x,d$y)

probs2 <- alpha_hat*d2(xo)/d1(xo)
hist(probs2)
summary(probs2)
identical(probs,probs2)
## turn probs and cl in ROC curve



plot(density(x2))
predict(density(x2),0:10)

sum(diff(out$fs[which.max(d2),]) > 0)
