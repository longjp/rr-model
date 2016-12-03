rm(list=ls())
library(Iso)
set.seed(1234)


## simulate two classes
n <- 1000
alpha <- 0.5
cl <- sample(1:2,n,replace=TRUE,prob=c(1-alpha,alpha))
x <- rnorm(n,mean=c(0,3)[cl],sd=1)
plot(density(x,bw="SJ"))


xs <- sort(x)

plot(xs,(1:n)/n)



cdf <- ((1:n)/n - (1-alpha)*pnorm(xs)) / alpha
out <- pava(cdf)
out[out<0] <- 0
out[out>1] <- 1

plot(xs,cdf,type='l',lwd=2)
points(xs,out,type='l',col='blue',lwd=2)
points(xs,pnorm(xs,mean=3),col="orange",lwd=2,type='l')

### run non-decreasing least squares

EstMixMdl <- function(data,gridsize=200){
    n <- length(data)
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
    ## Note: for Uniform(0,1) F_b(x) = x
    ## Usually would need to CHANGE this
    Fb <- pnorm(data.1)
    ## Compute the weights (= frequency/n) of the unique data values, i.e., dF_n
    Freq <- diff(c(0,Fn.1))
    distance <- rep(0,gridsize)
    distance[0]<- sqrt(t((Fn.1-Fb)^2)%*%Freq)
    for(i in 1:gridsize)
    {
        a <- i/gridsize
        ## Assumes a value of the mixing proportion
        F.hat <- (Fn.1-(1-a)*Fb)/a
        ## Computes the naive estimator of F_s
        F.is=pava(F.hat,Freq,decreasing=FALSE) ## Computes the Isotonic Estimator of F_s
        F.is[which(F.is<=0)]=0
        F.is[which(F.is>=1)]=1
        distance[i] <- a*sqrt(t((F.hat-F.is)^2)%*%Freq);
    }
    return(distance)
}

gridsize <- 200
alpha_grid <- (1:gridsize)/gridsize
d <- EstMixMdl(x,gridsize=gridsize)
d1 <- diff(c(d[1],d))
d2 <- diff(c(d1[1],d1))
d2 <- d2 / max(d2)
d <- d / max(d)
plot(alpha_grid,d,type='l',xaxs='i')
points(alpha_grid,d2,lty=2,type='l')
abline(v=alpha)
alpha_hat <- alpha_grid[which.max(d2)]
alpha_hat
alpha
