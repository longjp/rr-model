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

NewtonUpdate <- function(m,t,phi){
    gp <- SawPrime(t+phi)
    g <- Saw(t+phi)
    del <- -2*sum(m*gp) + 2*sum(g*gp)
    h <- 2*sum(gp*gp)
    return((phi - h^{-1}*del) %% 1)
}

n <- 100
e.sd <- 1
t <- runif(n)
phi <- runif(1)
m <- Saw(t+phi) + rnorm(n,e.sd)
plot(t,m)


N <- 1000
phis_grid <- (0:(N-1))/N
rss <- rep(0,N)
tm <- proc.time()
for(ii in 1:N){
    m.pred <- Saw(t+phis_grid[ii])
    rss[ii] <- sum((m - m.pred)^2)
}
proc.time() - tm
plot(phis_grid,rss)
abline(v=phi,lwd=2,col='grey')
abline(v=phis_grid[which.min(rss)],lwd=2)



phi.init <- runif(1)
NN <- 21
phis_new <- rep(0,NN)
phis_new[1] <- phi.init
tm <- proc.time()
for(ii in 2:NN){
    phis_new[ii] <- NewtonUpdate(m,t,phis_new[ii-1])
}
proc.time() - tm
abline(v=phis_new[1],col='red',lwd=2)
abline(v=phis_new[2:(NN-1)],col=alpha('red',.1))
abline(v=phis_new[NN],col='blue',lwd=2,lty=2)
## phis_new[NN] is close to phis_grid[which.min(rss)]

