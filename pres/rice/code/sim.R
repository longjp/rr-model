## makes plots for simulation results from paper
rm(list=ls())
source('sim_funcs.R')
library(ellipse)
library(xtable)
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
draw_sigma <- function(n){
    return(sample(c(.01,.1,1),n,replace=TRUE,prob=c(.05,.9,.05)))
}

## nonlinear component of f
g <- function(x) f(x) - beta1 - x*beta2

## construct response
draw_y <- function(x,sigma){
    y <- f(x) + rnorm(length(x),0,sigma)
}

## empirically compute (E[xx^T]^{-1}E[xf(x)])
x <- draw_x(100000)
X <- cbind(1,x)
ftemp <- matrix(f(x),nrow=length(x))
solve(t(X)%*%X,t(X)%*%ftemp)
## should approximately match beta1, beta2
beta1 <- -1/6
beta2 <- 1
## compute approximation to true delta
a1 <- mean(g(x)^2)
a2 <- mean(x*g(x)^2)
a4 <- mean(x^2*g(x)^2)
Ap <- matrix(c(a1,a2,a2,a4),nrow=2)
B <- solve(t(X)%*%X/length(x))
A <- B%*%Ap%*%B
Delta <- ComputeDelta(A,B) ## apply gamma function

## plot one realization of the data
cex.lab <- 2
cex.axis <- 1.6
x <- draw_x(n)
sigma <- draw_sigma(n)
y <- draw_y(x,sigma)
lmfit <- lm(y ~ x,weights=sigma^{-2})
pdf("../figs/data.pdf",height=6,width=8)
par(mar=c(5,5,1.5,1))
plot(0,0,xlim=c(0,1),ylim=c(min(beta1,0),max(beta1+beta2,1)),
     col=0,xlab="x",ylab="y",cex.lab=cex.lab,cex.axis=cex.axis,xaxs='i')
segments(x,y+sigma,x,y-sigma,col='grey')
abline(lmfit,lwd=3,col="blue")
abline(a=-1/6,b=1,lwd=3,col="orange")
points(x,y)
xpts <- (0:100)/100
points(xpts,xpts^2,type='l',col='black',lwd=3)
ex <- c(expression("f(x)=x"^2),expression(paste(x^T*beta*plain("=-1/6 + x"))),
        expression(paste(plain("WLS Fit with W=")*Sigma^-1)))
legend("topleft",ex,col=c("black","orange","blue"),lwd=3,cex=2)
dev.off()


## compare ols, wls, weighting with estimated delta, weighting with unknown variances
## store parameter estimates
wls <- matrix(0,nrow=N,ncol=2)
ols <- matrix(0,nrow=N,ncol=2)
dls <- matrix(0,nrow=N,ncol=2)
uls <- matrix(0,nrow=N,ncol=2)
deltahat <- rep(0,N)
## store asymptotic variance estimates
nuWLS1 <- array(0,dim=c(N,2,2))
nuOLS1 <- array(0,dim=c(N,2,2))
nuDLS1 <- array(0,dim=c(N,2,2))
nuWLS2 <- array(0,dim=c(N,2,2))
nuOLS2 <- array(0,dim=c(N,2,2))
nuDLS2 <- array(0,dim=c(N,2,2))
nuULS2 <- array(0,dim=c(N,2,2))
nuWLSOR <- array(0,dim=c(N,2,2))
nuOLSOR <- array(0,dim=c(N,2,2))
nuDLSOR <- array(0,dim=c(N,2,2))
print(paste0("beginning simulation"))
for(ii in 1:N){
    print(paste0("run ",ii," of ",N))
    ## draw data and compute wls and ols estimates
    x <- draw_x(n)
    sigma <- draw_sigma(n)
    y <- draw_y(x,sigma)
    ## compute ols and wls estimators
    wls.fit <- lm(y~x,weights=sigma^{-2})
    wls[ii,] <- wls.fit$coefficients
    ols.fit <- lm(y~x)
    ols[ii,] <- ols.fit$coefficients
    ## estimate B,A,delta and param estimates with adaptive weights
    Bhat <- ComputeBhat(x)
    Ahat <- ComputeAhat(y,x,sigma,Bhat)
    deltahat[ii] <- ComputeDelta(Ahat,Bhat)
    dls.fit <- lm(y~x,weights=(sigma^2 + deltahat[ii])^{-1})
    dls[ii,] <- dls.fit$coefficients
    ## estimates with unknown error variances
    m <- as.numeric(factor(sigma,labels=1:length(unique(sigma))))
    w_uls <- ComputeUnknownWeights(y,x,m,JJ=2)
    uls.fit <- lm(y~x,weights=w_uls)
    uls[ii,] <- uls.fit$coefficients
    ## compute asymptotic variance of all methods
    nuOLS1[ii,,] <- Nu1(Ahat,Bhat,sigma,rep(1,n),n)
    nuWLS1[ii,,] <- Nu1(Ahat,Bhat,sigma,sigma^{-2},n)
    nuDLS1[ii,,] <- Nu1(Ahat,Bhat,sigma,(sigma^2 + deltahat[ii])^{-1},n)
    nuOLSOR[ii,,] <- Nu1(A,B,sigma,rep(1,n),n)
    nuWLSOR[ii,,] <- Nu1(A,B,sigma,sigma^{-2},n)
    nuDLSOR[ii,,] <- Nu1(A,B,sigma,(sigma^2 + deltahat[ii])^{-1},n)
    nuOLS2[ii,,] <- Nu2(Bhat,ols.fit$residuals,rep(1,n),x)
    nuWLS2[ii,,] <- Nu2(Bhat,wls.fit$residuals,sigma^{-2},x)
    nuDLS2[ii,,] <- Nu2(Bhat,dls.fit$residuals,(sigma^2 + deltahat[ii])^{-1},x)
    nuULS2[ii,,] <- Nu2(Bhat,uls.fit$residuals,w_uls,x)
}


## these should all be 1, covariance matrices all psd
mean(apply(nuOLS1,1,function(x){min(eigen(x)$values) > 0}))
mean(apply(nuWLS1,1,function(x){min(eigen(x)$values) > 0}))
mean(apply(nuDLS1,1,function(x){min(eigen(x)$values) > 0}))
mean(apply(nuOLS2,1,function(x){min(eigen(x)$values) > 0}))
mean(apply(nuWLS2,1,function(x){min(eigen(x)$values) > 0}))
mean(apply(nuDLS2,1,function(x){min(eigen(x)$values) > 0}))
mean(apply(nuULS2,1,function(x){min(eigen(x)$values) > 0}))

## get fraction of times beta is in 95% CI
InCI <- function(theta,sigma,p,cip){
    c <- qchisq(cip,p)
    return((theta%*%solve(sigma)%*%t(theta) < c)[1,1])
}

InCIW <- function(thetas,sigmas,cip=0.95){
    if(length(dim(sigmas))==2){
        sigmas <- array(sigmas,dim=c(2,2,nrow(thetas)))
        sigmas <- aperm(sigmas,c(3,1,2))
    }
    thetas <- t(t(thetas) - c(beta1,beta2))
    p <- ncol(thetas)
    in_ci <- rep(FALSE,nrow(thetas))
    for(ii in 1:length(in_ci)){
        in_ci[ii] <- InCI(thetas[ii,,drop=FALSE],sigmas[ii,,],p,cip)
    }
    return(in_ci)
}

CItab <- matrix("-----",nrow=3,ncol=4)
CItab[1,1] <- mean(InCIW(wls,nuWLS1))
CItab[2,1] <- mean(InCIW(wls,nuWLS2))
CItab[3,1] <- mean(InCIW(wls,nuWLSOR))
CItab[1,2] <- mean(InCIW(ols,nuOLS1))
CItab[2,2] <- mean(InCIW(ols,nuOLS2))
CItab[3,2] <- mean(InCIW(ols,nuOLSOR))
CItab[1,3] <- mean(InCIW(dls,nuDLS1))
CItab[2,3] <- mean(InCIW(dls,nuDLS2))
CItab[3,3] <- mean(InCIW(dls,nuDLSOR))
CItab[2,4] <- mean(InCIW(uls,nuULS2))

rownames(CItab) <- c("$\\widehat{\\nu}_1$","$\\widehat{\\nu}_2$","$\\widehat{\\nu}_{OR}$")
colnames(CItab) <- c("WLS","OLS",
                     "$(\\Sigma + \\widehat{\\Delta})^{-1}$",
                     "$\\frac{\\Gamma(\\widehat{B})}{\\Gamma(\\widehat{B}^T\\widehat{C}_{m_i}\\widehat{B})}$")


cap <- "Fraction of times $\\beta$ is in 95\\% confidence region."
x.tab <- xtable(CItab,
                align="ccccc",
                caption=cap,
                label="tab:CI")
print(x.tab,
      file=paste0("../figs/CI.tex"),
      type='latex',
      include.rownames=TRUE,
      sanitize.text.function = function(x){x})

####### plot results
sig <- draw_sigma(10000) ## used for plotting asymptotic variance
cex.main <- 1.7
cex.lab <- 1.4
points_col <- "#00000020"
xlim <- range(c(ols[,1],wls[,1],dls[,1],uls[,1]))
ylim <- range(c(ols[,2],wls[,2],dls[,2],uls[,2]))

## ols
pdf("../figs/estimator_ols.pdf",width=6,height=5)
par(mar=c(5,5,1,1))
plot(ols,xlim=xlim,ylim=ylim,col=points_col,
     xlab=expression(hat(beta)[1]),
     ylab=expression(hat(beta)[2]),
     cex.lab=cex.lab)
points(ellipse(centre=c(beta1,beta2),Nu1(A,B,sig,rep(1,n),n),level=.95),
       type='l',col="red",lwd=2)
points(beta1,beta2,col="red",pch=19,cex=1.5)
dev.off()

## wls
pdf("../figs/estimator_wls.pdf",width=6,height=5)
par(mar=c(5,5,1,1))
plot(wls,xlim=xlim,ylim=ylim,col=points_col,
     xlab=expression(hat(beta)[1]),
     ylab=expression(hat(beta)[2]),
     cex.lab=cex.lab)
points(ellipse(centre=c(beta1,beta2),Nu1(A,B,sig,sig^{-2},n),level=.95),type='l',col="red",lwd=2)
points(beta1,beta2,col="red",pch=19,cex=1.5)
dev.off()

## adaptive, known variance
pdf("../figs/estimator_dls.pdf",width=6,height=5)
par(mar=c(5,5,1,1))
plot(dls,xlim=xlim,ylim=ylim,col=points_col,
     xlab=expression(hat(beta)[1]),
     ylab=expression(hat(beta)[2]),
     cex.lab=cex.lab)
points(ellipse(centre=c(beta1,beta2),Nu1(A,B,sig,(sig^2 + Delta)^{-1},n),level=.95),type='l',col="red",lwd=2)
points(beta1,beta2,col="red",pch=19,cex=1.5)
dev.off()

## adaptive, unknown variance
pdf("../figs/estimator_uls.pdf",width=6,height=5)
par(mar=c(5,5,1,1))
plot(uls,xlim=xlim,ylim=ylim,col=points_col,
     xlab=expression(hat(beta)[1]),
     ylab=expression(hat(beta)[2]),
     cex.lab=cex.lab)
points(ellipse(centre=c(beta1,beta2),Nu1(A,B,sig,(sig^2 + Delta)^{-1},n),level=.95),type='l',col="red",lwd=2)
points(beta1,beta2,col="red",pch=19,cex=1.5)
dev.off()



############# SAME AS ABOVE BUT BIGGER CAPTIONS FOR PRESENTATIONS
rmar <- 0.5
lmar <- 5.5
cex.lab <- 2
cex.axis <- 1.5
cex.main <- 3
xlim <- c(-.35,.05)
ylim <- c(.65,1.35)
height <- 6
width <- 6
## ols
pdf("../figs/pres_estimator_ols.pdf",width=width,height=height)
par(mar=c(5,lmar,5,rmar))
plot(ols,xlim=xlim,ylim=ylim,col=points_col,
     xlab=expression(hat(beta)[1]),
     ylab=expression(hat(beta)[2]),
     cex.lab=cex.lab,cex.axis=cex.axis,
     main="b) I",cex.main=cex.main)
points(ellipse(centre=c(beta1,beta2),Nu1(A,B,sig,rep(1,n),n),level=.95),
       type='l',col="red",lwd=2)
points(beta1,beta2,col="red",pch=19,cex=1.5)
dev.off()

## wls
pdf("../figs/pres_estimator_wls.pdf",width=width,height=height)
par(mar=c(5,lmar,5,rmar))
plot(wls,xlim=xlim,ylim=ylim,col=points_col,
     xlab=expression(hat(beta)[1]),
     ylab=expression(hat(beta)[2]),
     cex.lab=cex.lab,cex.axis=cex.axis,
     main=expression(plain("a) ")*Sigma^{-1}),cex.main=cex.main)
points(ellipse(centre=c(beta1,beta2),Nu1(A,B,sig,sig^{-2},n),level=.95),type='l',col="red",lwd=2)
points(beta1,beta2,col="red",pch=19,cex=1.5)
dev.off()

## adaptive, known variance
pdf("../figs/pres_estimator_dls.pdf",width=width,height=height)
par(mar=c(5,lmar,5,rmar))
plot(dls,xlim=xlim,ylim=ylim,col=points_col,
     xlab=expression(hat(beta)[1]),
     ylab=expression(hat(beta)[2]),
     cex.lab=cex.lab,cex.axis=cex.axis,
     main=expression(plain("c) ")*(Sigma + Gamma(hat(A))*Gamma(hat(B))^{-1})^{-1}),cex.main=cex.main)
points(ellipse(centre=c(beta1,beta2),Nu1(A,B,sig,(sig^2 + Delta)^{-1},n),level=.95),type='l',col="red",lwd=2)
points(beta1,beta2,col="red",pch=19,cex=1.5)
dev.off()

