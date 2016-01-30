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

## computes beta
ComputeBeta <- function(m,t,ampj,phij,phi,omega){
    X <- cbind(1,ampj*Saw(omega*t+phij+phi))
    B <- t(X)%*%X
    d <- t(X)%*%m
    z <- solve(B,d)
   return(z)
}

ComputeResiduals <- function(m,t,ampj,phij,phi,omega){
    X <- cbind(1,ampj*Saw(omega*t+phij+phi))
    z <- ComputeBeta(m,t,ampj,phij,phi,omega)
    return(m - X%*%z)
}

NewtonUpdate <- function(m,t,ampj,phij,params,omega){
    ## 1. condition on phi, omega, closed for update for amp,beta0
    beta0 <- params[1]
    amp <- params[2]
    phi <- params[3]
    beta <- ComputeBeta(m,t,ampj,phij,phi,omega)
    beta0 <- beta[1]
    amp <- beta[2]
    ## 2. condition on amp,beta0,omega, newton update phi (see sim.R for code)
    ## if amp = 0 we are at a (bad) stationary point, so choose random phase
    if(amp > 0){
        gp <- ampj*SawPrime(omega*t+phij+phi)
        g <- ampj*Saw(omega*t+phij+phi)
        del <- sum(gp*(g-(m-beta0)/amp))
        h <- sum(gp*gp)
        phi <- (phi - h^{-1}*del) %% 1
    } else {
        amp <- 0
        phi <- runif(1)
    }
    out <- c(beta0,amp,phi)
    return(out)
}

SawRss <- function(m,t,ampj,phij,omegas,NN=1){
    param <- c(0,1,runif(1))
    rss_min <- sum((m - mean(m))^2)
    rss <- rep(0,length(omegas))
    for(ii in 1:length(omegas)){
        for(jj in 1:NN){
            param <- NewtonUpdate(m,t,ampj,phij,param,omegas[ii])
        }
        X <- cbind(1,ampj*Saw(omegas[ii]*t+phij+param[3]))
        rss[ii] <- min(sum((m - X%*%param[1:2])^2),rss_min)
    }
    return(rss)
}
    
