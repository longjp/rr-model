### compare newton method to grid search
### for newton, condition on amplitude, beta0, update phase
### optimize only across phase

rm(list=ls())
library(scales)

## simulate from sawtooth with amp 1 and phase 0
##   t : evaluaation points
##  cc : fraction down
Saw <- function(t,cc=0.75){
    t <- t %% 1
    return(((-2/cc)*t + 1)*(t < cc) + ((2/(1-cc))*t - (1+cc)/(1-cc))*(t >= cc))
}

SawPrime <- function(t,cc=0.75){
    ## returns derivative of sawtooth evaluated at t
    t <- t %% 1
    return((-2/cc)*(t < cc) + ((2/(1-cc))*(t >= cc)))
}

## computes beta and projects onto amp > 0 axis
ComputeBeta <- function(m,t,phi,omega){
    X <- cbind(1,Saw(omega*t+phi))
    B <- t(X)%*%X
    d <- t(X)%*%m
    z <- solve(B,d)
    ## find best solution with amp > 0
    if(z[2] < 0){
        e <- c(0,1)
        q <- solve(B,e)
        z <- z - q*(sum(e*z)/sum(e*q))
        z[2] <- 0
    }
   return(z)
}

ComputeResiduals <- function(m,t,phi,omega){
    X <- cbind(1,Saw(omega*t+phi))
    z <- ComputeBeta(m,t,phi,omega)
    return(m - X%*%z)
}

NewtonUpdate <- function(m,t,params,omega){
    ## 1. condition on phi, omega, closed for update for amp,beta0
    beta0 <- params[1]
    amp <- params[2]
    phi <- params[3]
    beta <- ComputeBeta(m,t,phi,omega)
    beta0 <- beta[1]
    amp <- beta[2]
    ## 2. condition on amp,beta0,omega, newton update phi (see sim.R for code)
    ## if amp = 0 we are at a (bad) stationary point, so choose random phase
    if(amp > 0){
        gp <- SawPrime(omega*t+phi)
        g <- Saw(omega*t+phi)
        del <- sum(gp*(g-(m-beta0)/amp))
        h <- sum(gp*gp)
        phi <- (phi - h^{-1}*del) %% 1
    } else {
        phi <- runif(1)
    }
    out <- c(beta0,amp,phi)
    return(out)
}

omega_grid <- (1:1000)/1000 + 2

## construct data
n <- 25
e.sd <- 0.4
t <- runif(n)
phi <- runif(1)
amp <- rchisq(1,30)/30
beta0 <- rnorm(1)
omega <- sample(omega_grid,1)
m <- amp*Saw(omega*t+phi) + beta0 + rnorm(n,mean=0,sd=e.sd)
plot(t,m)

#### grid search to estimate omega and phase, closed for solution for beta0 and amp
N <- 500
phis_grid <- (0:(N-1))/N
tm <- proc.time()
rss <- rep(0,length(omega_grid))
for(ii in 1:length(rss)){
    print(ii)
    m.resid <- vapply(phis_grid,function(x){ComputeResiduals(m,t,x,omega_grid[ii])},rep(0,length(t))) ## faster not using lm
    rss[ii] <- min(colSums(m.resid^2))
}
proc.time() - tm
plot(omega_grid,rss)
abline(v=omega,col="black")

#### Newton optimization
rss <- rep(0,length(omega_grid))
NN <- 1
param <- c(0,1,runif(1))
tm <- proc.time()
for(ii in 1:length(omega_grid)){
    for(jj in 1:NN){
        param <- NewtonUpdate(m,t,param,omega_grid[ii])
    }
    X <- cbind(1,Saw(omega_grid[ii]*t+param[3]))
    rss[ii] <- sum((m - X%*%param[1:2])^2)
}
proc.time() - tm
points(omega_grid,rss,col='blue')
