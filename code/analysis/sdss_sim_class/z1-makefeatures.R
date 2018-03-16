### REWRITE UNCERTAINTY PIPELINE SO DO NOT COMPUTE FREQUENCY UNCERTAINTY

## compute uncertainties on all parameter estimates
## output sds (square root diagonal of hessian) along with features
rm(list=ls())

require(numDeriv)
load("../../fit_template/template_sdss.RData")
load("../../fit_template/template_des.RData")
source("../../fit_template/fit_template.R")
load("../../fit_template/template_sdss.RData")
source("../funcs.R")
source("../../common/funcs.R")

## load data
load("../../data/clean/sdss_sim_class.RData")
load("z0-fit.RData")


ComputeUncertainty <- function(coeffs,omega,lc,tem,use.errors=TRUE,use.dust=TRUE){
    if(use.dust){
        use.dust <- CheckNumberBands(lc)
    }
    CheckLC(lc)
    ## unpack coefficients
    tem <- CheckTemLC(tem,lc)
    dat <- AugmentData(lc,tem,use.errors)
    mag <- dat[[1]]$mag
    dust <- dat[[1]]$dust
    t <- dat[[1]]$time
    weights <- 1 / dat[[1]]$error^2
    nb <- dat[[2]]
    coeffs <- c(omega,coeffs)
    return(solve(hessian(ComputeRSS,coeffs,mag=mag,tem=tem,nb=nb,t=t,dust=dust,weights=weights)))
}

ComputeRSS <- function(coeffs,mag,tem,nb,t,dust,weights){
    omega <- coeffs[1]
    phi <- coeffs[5]
    coeffs <- coeffs[2:4]
    mag <- mag - rep.int(tem$abs_mag(1/omega,tem)[1,],nb)
    gammaf <- ConstructGamma(t,nb,phi,omega,tem$template_funcs)
    resid <- mag - coeffs[1] - coeffs[2]*dust - coeffs[3]*gammaf
    return(sum(weights*resid^2))
}



####### CHECK PERIODS REASONABLE
N <- sum(cl=="rr")
## fraction of times period is best
print("accuracies, top period:")
period_est <- period_est[,1] ## just use best fit period
print(paste("1%:",mean(abs((periods[1:N] - period_est[1:N])/periods[1:N]) < 0.01)))
print(paste("0.1%:",mean(abs((periods[1:N] - period_est[1:N])/periods[1:N]) < 0.001)))
print(paste("0.01%:",mean(abs((periods[1:N] - period_est[1:N])/periods[1:N]) < 0.0001)))

## fraction of times period is best on FULL
print("accuracies, top period:")
period_est_FULL <- period_est_FULL[,1] ## just use best fit period
print(paste("1%:",mean(abs((periods[1:N] - period_est_FULL[1:N])/periods[1:N]) < 0.01)))
print(paste("0.1%:",mean(abs((periods[1:N] - period_est_FULL[1:N])/periods[1:N]) < 0.001)))
print(paste("0.01%:",mean(abs((periods[1:N] - period_est_FULL[1:N])/periods[1:N]) < 0.0001)))


####### EXTRACT FEATURES
## a,ebv,phase,period,rss,variability measures,colors
phis <- (1:100)/100

## on full light curves, get features and uncertainties
coeffs_FULL <- matrix(0,nrow=length(tms_FULL),ncol=4)
for(ii in 1:length(tms_FULL)){
    tm <- tms_FULL[[ii]]
    lc <- TMtoLC(tm)
    p_est <- period_est_FULL[ii]
    omega <- 1 / p_est
    rssphi <- ComputeRSSPhase(lc,omega,tem_sdss,phis)
    phi <- phis[which.min(rssphi)]
    coeffs_FULL[ii,] <- ComputeCoeffsPhase(lc,omega,phi,tem_sdss)
}

coeffsu_FULL <- matrix(0,nrow=length(tms_FULL),ncol=5)
for(ii in 1:nrow(coeffsu_FULL)){
    u <- ComputeUncertainty(coeffs_FULL[ii,],1/period_est_FULL[ii],TMtoLC(tms_FULL[[ii]]),tem_sdss)
    coeffsu_FULL[ii,] <- diag(u)
}

## on poorly observed light curves, get features and uncertainties
coeffs <- matrix(0,nrow=length(tms),ncol=4)
for(ii in 1:length(tms)){
    tm <- tms[[ii]]
    lc <- TMtoLC(tm)
    p_est <- period_est[ii]
    omega <- 1 / p_est
    rssphi <- ComputeRSSPhase(lc,omega,tem_sdss,phis)
    phi <- phis[which.min(rssphi)]
    coeffs[ii,] <- ComputeCoeffsPhase(lc,omega,phi,tem_sdss)
}

coeffsu <- matrix(0,nrow=length(tms),ncol=4)
for(ii in 1:nrow(coeffsu)){
    u <- ComputeUncertainty(coeffs[ii,],1/period_est[ii],TMtoLC(tms[[ii]]),tem_sdss)
    coeffsu[ii,] <- diag(solve(solve(u)[2:5,2:5]))
}


sum(coeffsu<0)

## bind period estimates and coefficients
coeffs <- cbind(period_est,coeffs)
coeffs_FULL <- cbind(period_est_FULL[,1],coeffs_FULL)


## take square root of coeffsu
coeffsu <- sqrt(coeffsu)
coeffsu_FULL <- sqrt(coeffsu_FULL)



## QUESTIONS
## 1. negative uncertainty?
## 2. if amplitude = 0 (such as for non-rrl), what does uncertainty mean?

## compare LDA versus RF with original pipeline
## if competitive, implement irina's idea with this data
## if not, read up on kernalized lda, can you propagate feature uncertainty?
## read a few of irina's articles


## if we assume uncertainty is proportional to number of points, then easier from
## an implementation perspective


## if uncertainty is negative, input 1/4 range
r_over_4 <- apply(apply(coeffs,2,range),2,diff)/4
for(ii in 1:ncol(coeffsu)){
    coeffsu[is.na(coeffsu[,ii]),ii] <- r_over_4[ii]
}
r_over_4 <- apply(apply(coeffs_FULL,2,range),2,diff)/4
for(ii in 1:ncol(coeffsu_FULL)){
    coeffsu_FULL[is.na(coeffsu_FULL[,ii]),ii] <- r_over_4[ii]
}








## plot noisy versus clean features
plot(coeffs[,1],coeffs[,4],ylim=c(0,3))
segments(coeffs[,1] - coeffsu[,1],coeffs[,4],coeffs[,1] + coeffsu[,1],coeffs[,4],col="#00000020")
segments(coeffs[,1],coeffs[,4] - coeffsu[,4],coeffs[,1],coeffs[,4] + coeffsu[,4],col="#00000020")
dev.new()
plot(coeffs_FULL[,1],coeffs_FULL[,4],ylim=c(0,3))
segments(coeffs_FULL[,1] - coeffsu_FULL[,1],coeffs_FULL[,4],coeffs_FULL[,1] + coeffsu_FULL[,1],coeffs_FULL[,4],col="#00000020")
segments(coeffs_FULL[,1],coeffs_FULL[,4] - coeffsu_FULL[,4],coeffs_FULL[,1],coeffs_FULL[,4] + coeffsu_FULL[,4],col="#00000020")



temp <- rep("",length(npoints))
temp[npoints < median(npoints)] <- "poor"
clpoor <- paste0(cl,temp)



cols <- c("black","red","orange","blue")
names(cols) <- c("rr","not","rrpoor","notpoor")

pchs <- c(1:4)
names(pchs) <- c("rr","not","rrpoor","notpoor")



## plot noisy versus clean features
plot(coeffs[,1],coeffs[,4],ylim=c(0,3),col=cols[clpoor],pch=pchs[cl])
dev.new()
plot(coeffs_FULL[,1],coeffs_FULL[,4],ylim=c(0,3),col=cols[cl],pch=pchs[cl])


to_use <- npoints >= median(npoints)
plot(coeffs[to_use,1],coeffs[to_use,4],ylim=c(0,3),col=cols[cl[to_use]],pch=pchs[cl[to_use]])
dev.new()
to_use <- npoints >= median(npoints)
plot(coeffs[!to_use,1],coeffs[!to_use,4],ylim=c(0,3),col=cols[cl[!to_use]],pch=pchs[cl[!to_use]])






## 





## read in des data for same light curve
fname <- "LC_402316_des.dat"
lc_des <- read.table(fname,header=TRUE,stringsAsFactors=FALSE)
lc_des <- lc_des[,c(1,4,2,3)]
names(lc_des) <- c("time","band","mag","error")

## construct frequency grid
freq_min <- 1.0/0.95 ## max period of rrab is about 0.95 days 
freq_max <- 1.0/0.4 ## min period of rrab is about 0.4 days
max_phase_error <- .01 ## maximum phase error fraction, see documentation
freq_space <- (max_phase_error*4)/(max(lc_des[,1]) - min(lc_des[,1]))
omegas <- seq(from=freq_min,to=freq_max,by=freq_space)


## estimate period with des templates
rss_des <- FitTemplate(lc_des,omegas,tem_des)
omega_est_des <- omegas[which.min(rss_des)]
p_est_des <- 1/omega_est_des
coeffs_des <- ComputeCoeffs(lc_des,omega_est_des,tem_des)
names(coeffs_des) <- c("mu","ebv","amp","phase")


des_poor <- ComputeUncertainty(coeffs_des,omega_est_des,lc_des,tem_des)
des_poor
