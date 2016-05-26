## simulate from sawtooth with amp 1 and phase 0
##   t : evaluaation points
##  cc : fraction down
Saw <- function(t,cc){
    t <- t %% 1
    return(((-2/cc)*t + 1)*(t < cc) + ((2/(1-cc))*t - (1+cc)/(1-cc))*(t >= cc))
}

SawPrime <- function(t,cc){
    ## returns derivative of sawtooth evaluated at t
    t <- t %% 1
    return((-2/cc)*(t < cc) + ((2/(1-cc))*(t >= cc)))
}

## computes beta
ComputeBeta <- function(m,t,ampj,phij,phi,omega,cc){
    X <- cbind(1,ampj*Saw(omega*(t+phij)+phi,cc))
    B <- t(X)%*%X
    d <- t(X)%*%m
    z <- solve(B,d)
   return(z)
}

ComputeResiduals <- function(m,t,ampj,phij,phi,omega,cc){
    X <- cbind(1,ampj*Saw(omega*(t+phij)+phi,cc))
    z <- ComputeBeta(m,t,ampj,phij,phi,omega,cc)
    return(m - X%*%z)
}


GetParamsPhiGrid <- function(m,t,ampj,phij,omega,cc,N=500){
    phis_grid <- (0:(N-1))/N
    m.resid <- vapply(phis_grid,
                      function(x){ComputeResiduals(m,t,ampj,phij,x,omega,cc)},
                      rep(0,length(t)))
    phi <- phis_grid[which.min(colSums(m.resid^2))]
    betas <- ComputeBeta(m,t,ampj,phij,phi,omega,cc)
    params <- c(beta=betas[1],amp=betas[2],phi=phi)
    return(params)
}



NewtonUpdate <- function(m,t,ampj,phij,params,omega,cc){
    ## 1. condition on phi, omega, closed for update for amp,beta0
    beta0 <- params[1]
    amp <- params[2]
    phi <- params[3]
    beta <- ComputeBeta(m,t,ampj,phij,phi,omega,cc)
    beta0 <- beta[1]
    amp <- beta[2]
    ## 2. condition on amp,beta0,omega, newton update phi (see sim.R for code)
    ## if amp = 0 we are at a (bad) stationary point, so choose random phase
    if(amp > 0){
        gp <- ampj*SawPrime(omega*(t+phij)+phi,cc)
        g <- ampj*Saw(omega*(t+phij)+phi,cc)
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

SawRss <- function(m,t,ampj,phij,omegas,NN=1,cc=0.25){
    param <- c(0,1,runif(1))
    rss_min <- sum((m - mean(m))^2)
    rss <- rep(0,length(omegas))
    for(ii in 1:length(omegas)){
        for(jj in 1:NN){
            param <- NewtonUpdate(m,t,ampj,phij,param,omegas[ii],cc)
        }
        X <- cbind(1,ampj*Saw(omegas[ii]*(t+phij)+param[3],cc))
        rss[ii] <- min(sum((m - X%*%param[1:2])^2),rss_min)
    }
    return(rss)
}
    
