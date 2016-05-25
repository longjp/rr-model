### compare newton method to grid search
### for estimating phase, amplitude, and mean

rm(list=ls())
library(scales)

Saw <- function(t,cc=0.75){
    ## returns sawtooth evaluated at t
    t <- t %% 1
    return(((-2/cc)*t + 1)*(t < cc) + ((2/(1-cc))*t - (1+cc)/(1-cc))*(t >= cc))
}

SawPrime <- function(t,cc=0.75){
    ## returns derivative of sawtooth evaluated at t
    t <- t %% 1
    return((-2/cc)*(t < cc) + ((2/(1-cc))*(t >= cc)))
}

## ComputeRss <- function(m,t,phi,omega){
##     X <- cbind(1,Saw(omega*t+phi))
##     B <- t(X)%*%X
##     d <- t(X)%*%m
##     z <- solve(B,d)
##     ##z[2] <- max(z[2],0)
##     return(m - X%*%z)
## }

ComputeRss <- function(m,t,phi,omega){
    X <- cbind(1,Saw(omega*t+phi))
    B <- t(X)%*%X
    d <- t(X)%*%m
    z <- solve(B,d)
    if(z[2] < 0){
        e <- c(0,1)
        q <- solve(B,e)
        z <- z - q*(sum(e*z)/sum(e*q))
    }
    return(m - X%*%z)
}



NewtonUpdate <- function(m,t,params,omega){
    beta0 <- params[1]
    amp <- params[2]
    phi <- params[3]
    gp <- SawPrime(omega*t+phi)
    g <- Saw(omega*t+phi)
    X <- cbind(1,g)
    Hul <- 2*t(X)%*%X
    Hlr <- 2*amp^2*sum(gp*gp)
    Hoff <- c(2*amp*sum(gp),-2*sum((m-beta0)*gp)+4*amp*sum(g*gp))
    h <- rbind(cbind(Hul,Hoff),c(Hoff,Hlr))
    theta <- matrix(c(beta0,amp),nrow=2)
    dgdtheta <- 2*t(X)%*%X%*%theta - 2*t(X)%*%matrix(m,nrow=length(m))
    dgdphi <- -2*amp*sum((m-beta0)*gp) + 2*amp^2*sum(g*gp)
    del <- matrix(c(dgdtheta,dgdphi),nrow=3)
    params <- params - solve(h)%*%del
    params[3] <- params[3] %% 1
    params[2] <- max(0,params[2]) ## amplitude must be positive
    return(params)
}

## construct data
n <- 1000
e.sd <- 0.05
t <- runif(n)
phi <- runif(1)
amp <- rchisq(1,30)/30
beta <- rnorm(1)
omega <- 2.4
m <- amp*Saw(omega*t+phi) + beta + rnorm(n,mean=0,sd=e.sd)
plot(t,m)




## grid search to estimate omega and phase, closed for solution for beta0 and amp
N <- 1000
phis_grid <- (0:(N-1))/N
tm <- proc.time()
##m.resid <- vapply(phis_grid,function(x){lm(m~Saw(omega_grid[ii]*t+x))$residuals},rep(0,length(t)))
m.resid <- vapply(phis_grid,function(x){ComputeRss(m,t,x,omega)},rep(0,length(t))) ## faster not using lm
rss <- colSums(m.resid^2)
proc.time() - tm
plot(phis_grid,rss)
abline(v=phi,col="black")
abline(v=(phi+.5) %% 1,col="black")

#### Newton optimization
NN <- 5
phis <- matrix(0,nrow=NN,ncol=3)
phis[1,] <- c(0,1,runif(1))
tm <- proc.time()
for(jj in 1:(NN-1)){
    phis[jj+1,] <- NewtonUpdate(m,t,phis[jj,],omega)
}
proc.time() - tm
abline(v=phis[1,3],col="red")
abline(v=phis[NN,3],col="blue")


