rm(list=ls())
library(Iso)
set.seed(1234)


## simulate two classes
n <- 1000
alpha <- 0.5
cl <- sample(1:2,n,replace=TRUE,prob=c(1-alpha,alpha))
x <- rnorm(n,mean=c(0,3)[cl],sd=1)
plot(density(x,bw="SJ"))



xo <- sort(x)

plot(xo,(1:n)/n)



cdf <- ((1:n)/n - (1-alpha)*pnorm(xo)) / alpha
out <- pava(cdf)
out[out<0] <- 0
out[out>1] <- 1

plot(xo,cdf,type='l',lwd=2)
points(xo,out,type='l',col='blue',lwd=2)
points(xo,pnorm(xo,mean=3),col="orange",lwd=2,type='l')

### run non-decreasing least squares
EstMixMdl <- function(data,FbD,alpha_grid=(1:200)/200){
    ## Length of the data set
    data <- sort(data)
    ## Sorts the data set
    data.1 <- unique(data)
    ## Finds the unique data points
    Fn <- ecdf(data)
    ## Computes the empirical DF of the data
    Fn.1 <- Fn(data.1)
    ## Empirical DF of the data at the data points
    ## Calculate the known F_b at the data points
    Fb <- FbD(data.1)
    ## Compute the weights (= frequency/n) of the unique data values, i.e., dF_n
    Freq <- diff(c(0,Fn.1))
    distance <- rep(0,length(alpha_grid))
    F.isS <- matrix(0,nrow=length(alpha_grid),ncol=length(Fn.1))
    for(ii in 1:length(alpha_grid)){
        ## Assumes a value of the mixing proportion
        F.hat <- (Fn.1-(1-alpha_grid[ii])*Fb)/alpha_grid[ii]
        ## Computes the naive estimator of F_s
        F.is <- pava(F.hat,Freq,decreasing=FALSE) ## Computes the Isotonic Estimator of F_s
        F.is[F.is<=0] <- 0
        F.is[F.is>=1] <- 1
        F.isS[ii,] <- F.is
        distance[ii] <- alpha_grid[ii]*sqrt(sum(((F.hat-F.is)^2)*Freq))
    }
    return(list(distance=distance,alpha_grid=alpha_grid,fs=F.isS))
}

ComputeSecondDeriv <- function(x){
    return(c(0,0,diff(diff(x))))
}

EstimateParams <- function(data,FbD){
    data <- sort(data)
    out <- EstMixMdl(data,FbD)
    ix <- which.max(ComputeSecondDeriv(out$distance))
    alpha_hat <- out$alpha_grid[ix]
    Fs_est <- out$fs[ix,]
    xs <- data[diff(c(0,Fs_est)) > 0]
    fs <- density(xs)
    fs <- approxfun(fs$x,fs$y,rule=2)
    f <- density(data)
    f <- approxfun(f$x,f$y,rule=2)
    return(list(alpha_hat=alpha_hat,fs=fs,f=f))
}


FbD <- pnorm
params <- EstimateParams(x,FbD)
probs <- (params$alpha_hat*params$fs(x)) / params$f(x)
hist(probs)



ComputeROCCurve <- function(cl,probs){
    cl_pred <- cl[order(probs,decreasing=TRUE)]
    return(cbind(cumsum(-cl_pred + 1)/sum(cl==0),cumsum(cl_pred)/sum(cl==1)))
}

plot(ComputeROCCurve(cl-1,probs),type='l')


probs_true <- alpha*dnorm(x,mean=3)





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
