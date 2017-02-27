## makes plots for simulation results from paper
rm(list=ls())
library(ellipse)
set.seed(1234)

## constants
n <- 100 ## sample size
N <- 1000 ## number times to run sim

## functions for making data sets

## true relationship between y and x
f <- function(x) x^2

## distribution of predictors
draw_x <- function(n){
    return(runif(n))
}

## distribution of error variances
draw_sigma <- function(x){
    sigma <- rep(.1,length(x))
    sigma[x<.05] <- .01
    sigma[x>.95] <- 1
    return(sigma)
}

## nonlinear component of f
g <- function(x) f(x) - beta1 - x*beta2

## construct response
draw_y <- function(x,sigma){
    y <- f(x) + rnorm(length(x),0,sigma)
}

## empirically compute (E[xx^T]^{-1}E[xf(x)])
x <- draw_x(100000)
sig <- draw_sigma(x)
X <- cbind(1,x)
ftemp <- matrix(f(x),nrow=length(x))
solve(t(X)%*%X,t(X)%*%ftemp)
## should approximately match beta1, beta2
beta1 <- -1/6
beta2 <- 1
## compute (very good) approximation to sandwich covariance
a1 <- mean(g(x)^2 + sig^2)
a2 <- mean(x*(g(x)^2 + sig^2))
a4 <- mean(x^2*(g(x)^2 + sig^2))
meat <- matrix(c(a1,a2,a2,a4),nrow=2)
bread <- solve(t(X)%*%X/nrow(X))
sand <- bread%*%meat%*%bread
## compute what wls is converging to
a1 <- mean(sig^{-2})
a2 <- mean(x*sig^{-2})
a4 <- mean(x^2*sig^{-2})
T1 <- matrix(c(a1,a2,a2,a4),nrow=2)
T2 <- matrix(colMeans((X*sig^{-2})*f(x)),nrow=2)
beta_wls <- solve(T1)%*%T2


## compare ols, wls, weighting with estimated delta, weighting with unknown variances
## store parameter estimates
wls <- matrix(0,nrow=N,ncol=2)
ols <- matrix(0,nrow=N,ncol=2)
print(paste0("beginning simulation"))
for(ii in 1:N){
    print(paste0("run ",ii," of ",N))
    ## draw data and compute wls and ols estimates
    x <- draw_x(n)
    sigma <- draw_sigma(x)
    y <- draw_y(x,sigma)
    ## compute ols and wls estimators
    wls.fit <- lm(y~x,weights=sigma^{-2})
    wls[ii,] <- wls.fit$coefficients
    ols.fit <- lm(y~x)
    ols[ii,] <- ols.fit$coefficients
}

####### plot results
cex.main <- 1.7
cex.lab <- 1.4
points_col <- "#00000020"
xlim <- range(c(ols[,1],wls[,1]))
ylim <- range(c(ols[,2],wls[,2]))

## ols
pdf("figs/estimator_ols_dep.pdf",width=6,height=5)
par(mar=c(5,5,1,1))
plot(ols,xlim=xlim,ylim=ylim,col=points_col,
     xlab=expression(hat(beta)[1]),
     ylab=expression(hat(beta)[2]),
     cex.lab=cex.lab)
points(ellipse(centre=c(beta1,beta2),sand/n,level=.95),
       type='l',col="red",lwd=2)
points(beta1,beta2,col="red",pch=19,cex=1.5)
dev.off()

## wls
pdf("figs/estimator_wls_dep.pdf",width=6,height=5)
par(mar=c(5,5,1,1))
plot(wls,xlim=xlim,ylim=ylim,col=points_col,
     xlab=expression(hat(beta)[1]),
     ylab=expression(hat(beta)[2]),
     cex.lab=cex.lab)
points(beta1,beta2,col="red",pch=19,cex=1.5)
points(beta_wls[1],beta_wls[2],col="orange",pch=4,cex=1.5,lwd=3)
dev.off()
