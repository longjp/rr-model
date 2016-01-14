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

NewtonUpdate <- function(m,t,phi,amp,beta0,omega){
    gp <- SawPrime(omega*t+phi)
    g <- Saw(omega*t+phi)
    X <- cbind(1,g)
    Hul <- 2*t(X)%*%X
    Hlr <- 2*amp^2*sum(g*gp)
    Hoff <- c(2*amp*sum(g),-2*sum((m-beta0)*gp)+4*amp*sum(g*gp))
    h <- rbind(cbind(Hul,Hoff),c(Hoff,Hlr))
    Beta <- matrix(c(beta0,amp),nrow=2)
    dgdBeta <- 2*t(X)%*%X%*%Beta - 2*t(X)%*%m
    dgdphi <- -2*amp*sum((m-beta0)*gp) + 2*amp^2*sum(g*gp)
    del <- c(dgdBeta,dgdphi)
    theta <- c(beta0,amp,phi) - solve(h)*del
    theta[3] <- theta[3] %% 1
    return(theta)
}

omega_grid <- (1:1000)/1000 + 2

## construct data
n <- 1000
e.sd <- 0.05
t <- runif(n)
phi <- runif(1)
amp <- rchisq(1,30)/30
beta <- rnorm(1)
omega <- sample(omega_grid,1)
m <- amp*Saw(omega*t+phi) + beta + rnorm(n,mean=0,sd=e.sd)
plot(t,m)



compute_rss <- function(m,t,phi,omega){
    X <- cbind(1,Saw(omega*t+phi))
    B <- t(X)%*%X
    d <- t(X)%*%m
    z <- solve(B,d)
    return(m - X%*%z)
}



## grid search to estimate omega and phase
N <- 100
phis_grid <- (0:(N-1))/N
rss <- rep(0,length(omega_grid))
tm <- proc.time()
for(ii in 1:length(omega_grid)){
    ##m.resid <- vapply(phis_grid,function(x){lm(m~Saw(omega_grid[ii]*t+x))$residuals},rep(0,length(t)))
    m.resid <- vapply(phis_grid,function(x){compute_rss(m,t,x,omega_grid[ii])},rep(0,length(t))) ## faster not using lm
    rss[ii] <- min(colSums(m.resid^2))
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


