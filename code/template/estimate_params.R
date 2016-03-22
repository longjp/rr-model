rm(list=ls())
source("rrab_fit.R")
load("tms_params.RData")
load("make_template.RData")

temp_time <- seq(0,1,length.out=ncol(tem$templates))





ii <- 1
dat <- AugmentData(tms[[ii]],tem$dust,tem$betas)


par(mfcol=c(3,1))
plot((dat[[1]]$time %% param$period[ii]) / param$period[ii],dat[[1]]$mag - dat[[1]]$dust*param$d[ii] - param$alpha[ii])
plot((dat[[1]]$time %% param$period[ii]) / param$period[ii],dat[[1]]$mag - dat[[1]]$dust*param$d[ii])
plot((tms[[ii]]$time %% param$period[ii]) / param$period[ii],tms[[ii]]$mag)




m <- dat[[1]]$mag
dust <- dat[[1]]$dust
t <- dat[[1]]$time
nb <- dat[[2]]

phi <- param$phase[ii]
per <- param$period[ii]
omega <- 1/per

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




gammaf <- ConstructGamma(t,nb,phi,omega,template_funcs)
pred <- gammaf*param$a[ii] + param$d[ii]*dat[[1]]$dust + rep(param$alpha[ii],length(gammaf)) + rep.int(tem$betas,nb)

tms[[ii]] <- tms[[ii]][order(tms[[ii]][,2]),]

plot((tms[[ii]][,1] %% per) / per,tms[[ii]][,3])
points((tms[[ii]][,1] %% per) / per,pred,col='red')







## simple test of newton - given correct period can
## it find other parameters - answer yes
alpha <- param$alpha[ii]
d <- param$d[ii]
a <- param$a[ii]
param$phase[ii]
phi <- runif(1)
phi
for(jj in 1:5){
    out <- NewtonUpdate(phi,omega,m,t,dust,nb,template_funcs,templated_funcs)
    print(out)
    phi <- out[4]
}


mp <- m - alpha - d*dust

GammaError <- function(phi){
    gammaf <- ConstructGamma(t,nb,phi,omega,template_funcs)
    return(sum((mp - a*gammaf)^2))
}

phi_grid <- seq(0,1,length.out=100)
phi_grid_error <- vapply(phi_grid,GammaError,c(0))
plot(phi_grid,phi_grid_error)
abline(v=param$phase[ii])
