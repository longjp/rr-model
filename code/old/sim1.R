### compare newton method to grid search
### for estimating phase

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

NewtonUpdate <- function(m,t,phi,omega){
    gp <- SawPrime(omega*t+phi)
    g <- Saw(omega*t+phi)
    del <- sum(gp*(g-m))
    h <- sum(gp*gp)
    ## del <- -2*sum(m*gp) + 2*sum(g*gp)
    ## h <- 2*sum(gp*gp)
    return((phi - h^{-1}*del) %% 1)
}

omega_grid <- (1:1000)/1000 + 2

## construct data
n <- 20
e.sd <- 0.25
t <- runif(n)
phi <- runif(1)
omega <- sample(omega_grid,1)
m <- Saw(omega*t+phi) + rnorm(n,mean=0,sd=e.sd)
plot(t,m)

## grid search to estimate omega and phase
N <- 100
phis_grid <- (0:(N-1))/N
rss <- rep(0,length(omega_grid))
tm <- proc.time()
for(ii in 1:length(omega_grid)){
    m.pred <- vapply(phis_grid,function(x){Saw(omega_grid[ii]*t+x)},rep(0,length(t)))
    rss[ii] <- min(colSums((m-m.pred)^2))
}
proc.time() - tm
plot(omega_grid,rss)
abline(v=omega,col="black")

#### Newton optimization
rss <- rep(0,length(omega_grid))
NN <- 1 ## number of Newton steps for each updated frequency
phi_temp <- runif(1) ## random initial phase
tm <- proc.time()
for(ii in 1:length(omega_grid)){
    for(jj in 1:NN){
        phi_temp <- NewtonUpdate(m,t,phi_temp,omega_grid[ii])
    }
    rss[ii] <- sum((Saw(omega_grid[ii]*t+phi_temp) - m)^2)
}
proc.time() - tm
points(omega_grid,rss,col='red')


