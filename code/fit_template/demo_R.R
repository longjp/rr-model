
### load template fitting code and templates
source("fit_template.R") ## fits model
load("template_sdss.RData") ## template for sdss filters, in tem_sdss
load("template_des.RData") ## template for des filters, in tem_des
ls()

## read in sdss light curve and plot
fname <- "LC_4099.dat"
lc <- read.table(fname,stringsAsFactors=FALSE)
names(lc) <- c("time","band","mag","error")



## colors for plotting
options(repr.plot.width=6,repr.plot.height=4)
bandpch <- 1:6
names(bandpch) <- c("u","g","r","i","z","Y")
bandcol <- c("dodgerblue3","green","red",
             "mediumorchid1","black","peachpuff4")
names(bandcol) <- c("u","g","r","i","z","Y")

## plot raw light curve
par(mar=c(5,5,1,1))
plot(lc$time,lc$mag,col=bandcol[lc$band],pch=bandpch[lc$band],
     ylim=rev(range(lc$mag)),
     xlab="time",ylab="magnitude")
segments(lc$time,lc$mag+lc$error,lc$time,lc$mag-lc$error,col='grey')


## plots template fit on light curve, used later
plotLC <- function(lc,p_est,coeffs,tem){
    lc1 <- lc
    lc1[,1] <- (lc$time %% p_est)/p_est
    lc2 <- lc1
    lc2[,1] <- lc1[,1] + 1
    lc_temp <-rbind(lc1,lc2)
    plot(lc_temp$time,lc_temp$mag,
         col=bandcol[lc_temp$band],pch=bandpch[lc_temp$band],
         ylim=rev(range(lc_temp$mag)),
         xlab="phase",ylab="magnitude",
         xlim=c(0,2),xaxs='i')
    segments(lc_temp$time,
             lc_temp$mag+lc_temp$error,
             lc_temp$time,
             lc_temp$mag-lc_temp$error,col='grey')
    ti <- seq(0,p_est,length.out=100)
    ti <- c(ti,ti+p_est)
    m <- PredictAllBand(ti,1/p_est,coeffs,tem)
    for(ii in 1:ncol(m)){
        points(ti/p_est,m[,ii],type='l',col=bandcol[colnames(m)[ii]])
    }
}


## construct frequency grid
freq_min <- 1.0/0.95 ## max period of rrab is about 0.95 days 
freq_max <- 1.0/0.4 ## min period of rrab is about 0.4 days
max_phase_error <- .02 ## maximum phase error fraction, see documentation
freq_space <- (max_phase_error*4)/(max(lc[,1]) - min(lc[,1]))
omegas <- seq(from=freq_min,to=freq_max,by=freq_space)
length(omegas) ## longer omegas, more computation time, but less risk of missing rss min

## compute model fit at each frequency (low rss = good fit)
rss <- FitTemplate(lc,omegas,tem_sdss)
omega_est <- omegas[which.min(rss)]

## obtain parameter estimates
p_est <- 1/omega_est
p_est
coeffs <- ComputeCoeffs(lc,omega_est,tem_sdss) ## at fixed frequency, determines mu, ebv, amp, phase
names(coeffs) <- c("mu","ebv","amp","phase")
coeffs
## correct values near:
## > coeffs
##          mu         ebv         amp       phase 
## 16.17423839  0.05122597  0.58949080  0.24304752 
## p_est
## [1] 0.6417558

## check rss as function of period, hopefully single, unambiguous global min
plot(1/omegas,rss,xlab="period",ylab="rss")
abline(v=p_est)



coeffs <- c(omega_est,coeffs)
names(coeffs)[1] <- "omega"


library(numDeriv)

sdss_well <- solve(hessian(ComputeRSS,coeffs,lc=lc,tem=tem_sdss))
sqrt(diag(sdss_well))







## plot folded light curve with best fit
par(mar=c(5,5,1,1))
plotLC(lc,p_est,coeffs,tem_sdss)

## check that phase determined by NewtonUpdate is close to best phase
## by performing a grid search on phase, comparing min of grid
## search to min found by newton
phis <- (1:100)/100
rss_phi <- ComputeRSSPhase(lc,omega_est,tem_sdss,phis=phis)
par(mar=c(5,5,1,1))
plot(phis,rss_phi)
abline(v=coeffs[4])


ebv <- coeffs[2]
## remove dust
bands <- names(tem_sdss$dust)
for(ii in 1:length(bands)){
    lc$mag[lc$band==bands[ii]] <- lc$mag[lc$band==bands[ii]] - tem_sdss$dust[ii]*ebv
}

## compute RSS, get parameters
rss <- FitTemplate(lc,omegas,tem_sdss,use.dust=FALSE)
omega_est <- omegas[which.min(rss)]
p_est <- 1/omega_est
coeffs <- ComputeCoeffs(lc,omega_est,tem_sdss,use.dust=FALSE)
names(coeffs) <- c("mu","ebv","amp","phase")

## view rss
par(mar=c(5,5,1,1))
plot(1/omegas,rss,xlab="period",ylab="rss")
abline(v=p_est)

## plot folded light curve with best fit
par(mar=c(5,5,1,1))
plotLC(lc,p_est,coeffs,tem_sdss)

## check that phase determine by NewtonUpdate is close to best phase
## by performing a grid search
phis <- (1:100)/100
rss_phi <- ComputeRSSPhase(lc,omega_est,tem_sdss,phis=phis,use.dust=FALSE)
par(mar=c(5,5,1,1))
plot(phis,rss_phi)
abline(v=coeffs[4])

## add lots of error to selected observations
ix <- sample(1:nrow(lc),floor(nrow(lc)*.6))
lc[ix,3] <- lc[ix,3] + rnorm(length(ix),mean=0,sd=3)
lc[ix,4] <- sqrt(lc[ix,4]^2 + 3^2)

## plot raw light curve
colpch <- 1:5
names(colpch) <- unique(lc$band)
par(mar=c(5,5,1,1))
plot(lc$time,lc$mag,col=colpch[lc$band],pch=colpch[lc$band],
     ylim=rev(range(lc$mag)),
     xlab="time",ylab="magnitude")
segments(lc$time,lc$mag+lc$error,lc$time,lc$mag-lc$error,col='grey',lwd=0.5)

###### a) use.errors=TRUE
rss <- FitTemplate(lc,omegas,tem_sdss)
omega_est <- omegas[which.min(rss)]
p_est <- 1/omega_est
coeffs <- ComputeCoeffs(lc,omega_est,tem_sdss)
names(coeffs) <- c("mu","ebv","amp","phase")
coeffs
p_est

## view rss, still clear min
par(mar=c(5,5,1,1))
plot(1/omegas,rss,xlab="period",ylab="rss")
abline(v=p_est)

###### b) use.errors=FALSE, should get poor parameter estimates
rss <- FitTemplate(lc,omegas,tem_sdss,use.errors=FALSE)
omega_est <- omegas[which.min(rss)]
p_est <- 1/omega_est
coeffs <- ComputeCoeffs(lc,omega_est,tem_sdss,use.errors=FALSE)
names(coeffs) <- c("mu","ebv","amp","phase")
coeffs
p_est



## view rss, very messy
par(mar=c(5,5,1,1))
plot(1/omegas,rss,xlab="period",ylab="rss")
abline(v=p_est)






## read in sdss data
fname <- "LC_402316.dat"
lc_sdss <- read.table(fname,stringsAsFactors=FALSE)
names(lc_sdss) <- c("time","band","mag","error")

## fit model
rss_sdss <- FitTemplate(lc_sdss,omegas,tem_sdss)
omega_est_sdss <- omegas[which.min(rss_sdss)]
p_est_sdss <- 1/omega_est_sdss
coeffs_sdss <- ComputeCoeffs(lc_sdss,omega_est_sdss,tem_sdss)
names(coeffs_sdss) <- c("mu","ebv","amp","phase")
coeffs_sdss

## plot folded light curve with best fit
par(mar=c(5,5,1,1))
plotLC(lc_sdss,p_est_sdss,coeffs_sdss,tem_sdss)


## plot RSS as function of phase
phis <- (1:100)/100
rssphase <- ComputeRSSPhase(lc_sdss,omega_est_sdss,tem_sdss,phis=phis)
plot(phis,rssphase)
abline(v=coeffs_sdss[4])

## compute hessian
## coeffs_sdss <- c(omega_est_sdss,coeffs_sdss)
## names(coeffs_sdss)[1] <- "omega"
coeffs_sdss <- c(omega_est_sdss,coeffs_sdss)
names(coeffs_sdss)[1] <- "omega"
coeffs_sdss
sdss_well <- solve(hessian(ComputeRSS,coeffs_sdss,lc=lc_sdss,tem=tem_sdss))
sdss_well
sqrt(diag(sdss_well))


## computes the RSS at fit
## since the RSS is the log likelihood finding the Hessian of 
## coeffs = [omega,mu,ebv,amp,phase]
## ComputeRSS <- function(coeffs,omega,phi,lc,tem,use.errors=TRUE,use.dust=TRUE){
##     if(use.dust){
##         use.dust <- CheckNumberBands(lc)
##     }
##     CheckLC(lc)
##     ## unpack coefficients
##     tem <- CheckTemLC(tem,lc)
##     dat <- AugmentData(lc,tem,use.errors)
##     m <- dat[[1]]$mag
##     dust <- dat[[1]]$dust
##     t <- dat[[1]]$time
##     weights <- 1 / dat[[1]]$error^2
##     nb <- dat[[2]]
##     m <- m - rep.int(tem$abs_mag(1/omega,tem)[1,],nb)
##     rss_max <- sum((lm(m~dust,weights=weights)$residuals^2)*weights)
##     gammaf <- ConstructGamma(t,nb,phi,omega,tem$template_funcs)
##     resid <- m - coeffs[1] - coeffs[2]*dust - coeffs[3]*gammaf
##     ##return(min(sum(weights*resid^2),rss_max))
##     return(sum(weights*resid^2))
## }



######## WORKS FROM HERE
### TODO: 1) more testing of method
###       2) use / do not use dust
require(numDeriv)
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







## read in des data for same light curve
fname <- "LC_402316_des.dat"
lc_des <- read.table(fname,header=TRUE,stringsAsFactors=FALSE)
lc_des <- lc_des[,c(1,4,2,3)]
names(lc_des) <- c("time","band","mag","error")

## estimate period with des templates
rss_des <- FitTemplate(lc_des,omegas,tem_des)
omega_est_des <- omegas[which.min(rss_des)]
p_est_des <- 1/omega_est_des
coeffs_des <- ComputeCoeffs(lc_des,omega_est_des,tem_des)
names(coeffs_des) <- c("mu","ebv","amp","phase")


des_poor <- ComputeUncertainty(coeffs_des,omega_est_des,lc_des,tem_des)
des_poor




library(microbenchmark)
?microbenchmark
out <- microbenchmark(hessian(ComputeRSS,coeffs_des,lc=lc_des,tem=tem_des))
print(out)


x <- rnorm(1000)
out <- microbenchmark(sum(x))
summary(out)

## plot folded light curve with best fit
par(mar=c(5,5,1,1))
plotLC(lc_des,p_est_des,coeffs_des,tem_des)

## plot rss functions for sdss and des
## DES minima still clear
par(mfcol=c(2,1),mar=c(5,5,1,1))
plot(1/omegas,rss_sdss,xlab="period",ylab="rss sloan data")
abline(v=p_est_sdss)
plot(1/omegas,rss_des,xlab="period",ylab="rss des data")
abline(v=p_est_des)


p_est_des
p_est_sdss


