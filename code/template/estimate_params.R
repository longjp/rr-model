rm(list=ls())
source("rrab_fit.R")
load("tms_params.RData")
load("make_template.RData")

temp_time <- seq(0,1,length.out=ncol(tem$templates))





ii <- 1
dat <- AugmentData(tms[[ii]],tem$dust,tem$betas)


par(mfcol=c(2,1))
plot((dat[[1]]$time %% param$period[ii]) / param$period[ii],dat[[1]]$mag - dat[[1]]$dust*param$d[ii])
plot((tms[[ii]]$time %% param$period[ii]) / param$period[ii],tms[[1]]$mag)



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
        amp <- 0
        phi <- runif(1)
    }
    out <- c(beta0,amp,phi)
    return(out)
}



t <- dat[[1]]$time
nb <- dat[[2]]
phi <- param$phase[ii]
per <- param$period[ii]
m <- dat[[1]]$mag
dust <- dat[[1]]$dust

gammaf <- ConstructGamma(t,nb,phi,per,tem$templates,temp_time)


pred <- gammaf*param$a[ii] + param$d[ii]*dat[[1]]$dust + rep(param$alpha[ii],length(gammaf)) + rep.int(tem$betas,nb)

tms[[ii]] <- tms[[ii]][order(tms[[ii]][,2]),]

plot((tms[[ii]][,1] %% per) / per,tms[[ii]][,3])
points((tms[[ii]][,1] %% per) / per,pred,col='red')


templates <- tem$templates
templatesd <- tem$templatesd


template_funcs <- list()
for(jj in 1:nrow(templates)){
    template_funcs[[jj]] <- approxfun(temp_time,templates[jj,])
}
templated_funcs <- list()
for(jj in 1:nrow(templatesd)){
    templated_funcs[[jj]] <- approxfun(temp_time,templatesd[jj,])
}





gammaf <- ConstructGamma(t,nb,phi,per,template_funcs)
est <- ComputeBeta(m,dust,gammaf)
alpha <- est["alpha"]
a <- est["a"]
d <- est["d"]
if(a > 0){
    ## write newton update here
    gp <- SawPrime(omega*t+phi)
    g <- Saw(omega*t+phi)
    del <- sum(gp*(g-(m-beta0)/amp))
    h <- sum(gp*gp)
    phi <- (phi - h^{-1}*del) %% 1
} else {
    a <- 0
    phi <- runif(1)
}
out <- c(alpha,d,a,phi)


