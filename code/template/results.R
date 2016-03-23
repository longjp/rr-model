rm(list=ls())
load("tms_params.RData")
load("make_template.RData")
load("estimate_params.RData")
source("rrab_fit.R")

N <- 10

print("accuracies:")
print(paste("1%:",mean(abs((param$period[1:N] - period_est)/param$period[1:N]) < 0.01)))
print(paste("0.1%:",mean(abs((param$period[1:N] - period_est)/param$period[1:N]) < 0.001)))
print(paste("0.01%:",mean(abs((param$period[1:N] - period_est)/param$period[1:N]) < 0.0001)))


temp_time <- seq(0,1,length.out=ncol(tem$templates))
tem$template_funcs <- list()
for(jj in 1:nrow(tem$templates)){
    tem$template_funcs[[jj]] <- approxfun(temp_time,tem$templates[jj,])
}
tem$templated_funcs <- list()
for(jj in 1:nrow(tem$templatesd)){
    tem$templated_funcs[[jj]] <- approxfun(temp_time,tem$templatesd[jj,])
}



ComputeCoeffs <- function(tm,omega,tem,NN=10){
    dat <- AugmentData(tm,tem$dust,tem$betas)
    m <- dat[[1]]$mag
    dust <- dat[[1]]$dust
    t <- dat[[1]]$time
    nb <- dat[[2]]

    coeffs <- c(0,0,0,runif(1))
    for(jj in 1:NN){
        coeffs <- NewtonUpdate(coeffs[4],omega,m,t,dust,nb,tem$template_funcs,tem$templated_funcs)
    }
    return(coeffs)
}
 

## plot all bands with best fit parameters
##for(ii in 1:length(p_est)){
ii <- 1
tm <- tms[[ii]]
omega <- 1/period_est[ii]
coeffs <- ComputeCoeffs(tm,omega,tem)


##pdf(paste0("figs/",ii,".pdf"),height=12,width=8)
dat <- AugmentData(tm,tem$dust,tem$betas)



######## FIGURE OUT HOW TO PLOT ALL LIGHT CURVES

band_names <- unique(tm$band)
par(mar=c(5,3,3,1),mfcol=c(5,1))
jj <- 1
##for(jj in 1:length(band_names)){
tm_temp <- dat[[1]][1:nb
t <- ((1:N)/N)/omega
ampj <- tm_temp$ampj[1]
phij <- tm_temp$phij[1]
pred <- params["beta"] + params["amp"]*ampj*Saw(omega*(t+phij)+params["phi"],cc)
plot(t * omega,pred,type='l',xlab="Phase",ylab="Mag",main=band_names[jj])
points((tm_temp$t %% (1/omega))*omega,tm_temp$m)
##}
##dev.off()
##}
