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





m <- dat[[1]]$mag
dust <- dat[[1]]$dust
t <- dat[[1]]$time

nb <- dat[[2]]

phi <- param$phase[ii]
per <- param$period[ii]

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



## work on this
NewtonUpdate <- function(m,t,dust,nb,omega,template_funcs,templated_funcs){
    
    gammaf <- ConstructGamma(t,nb,phi,omega,template_funcs)
    est <- ComputeBeta(m,dust,gammaf)
    alpha <- est["alpha"]
    a <- est["a"]
    d <- est["d"]
    if(a > 0){
        ## write newton update here
        gammafd <- ConstructGamma(t,nb,phi,per,templated_funcs)
        mp <- m - alpha - d*dust
        del <- sum(gammafd*(mp-a*gammaf))
        h <- a*sum(gp*gp)
        phi <- (phi + h^{-1}*del) %% 1
    } else {
        a <- 0
        phi <- runif(1)
    }
    out <- c(alpha,d,a,phi)


