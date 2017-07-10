### demonstrates how to use some functions in fit_template.R
### PART 1 below is the most useful, other parts are more code checking
rm(list=ls())
source("fit_template.R")
load("template.RData")

## read in data and plot
fname <- "LC_4099.dat"
lc <- read.table(fname,stringsAsFactors=FALSE)
names(lc) <- c("time","band","mag","error")

## colors for plotting
bandpch <- 1:6
names(bandpch) <- c("u","g","r","i","z","Y")
bandcol <- c("dodgerblue3","green","red",
             "mediumorchid1","black","peachpuff4")
names(bandcol) <- c("u","g","r","i","z","Y")

## plot raw light curve
plot(lc$time,lc$mag,col=bandcol[lc$band],pch=bandpch[lc$band],
     ylim=rev(range(lc$mag)),
     xlab="time",ylab="magnitude")
segments(lc$time,lc$mag+lc$error,lc$time,lc$mag-lc$error,col='grey')

## true period of source is 0.6417558

## makes nice plot
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
    ti <- (1:100)/100
    ti <- c(ti,ti+1)
    m <- PredictAllBand(ti,1,coeffs,tem)
    for(ii in 1:ncol(m)){
        points(ti,m[,ii],type='l',col=bandcol[colnames(m)[ii]])
    }
}

##
## PART 1: test if model is working in the most standard mode
##         (fitting for dust and using photometric errors)
##         


## fit template model and obtain coefficients
omegas <- seq(from=1.0,to=5.0,by=0.1/4000.0)
rss <- FitTemplate(lc,omegas,tem)
omega_est <- omegas[which.min(rss)]
p_est <- 1/omega_est
coeffs <- ComputeCoeffs(lc,omega_est,tem)
names(coeffs) <- c("mu","ebv","amp","phase")
## correct values near:
## coeffs
##[1] 16.0914314  0.1071064  0.5444982  0.3935931
## view rss

plot(1/omegas,rss,xlab="period",ylab="rss")
abline(v=p_est)

## plot folded light curve with best fit
plotLC(lc,p_est,coeffs,tem)

## check that phase determined by NewtonUpdate is close to best phase
## by performing a grid search on phase, comparing min of grid
## search to min found by newton
phis <- (1:100)/100
rss_phi <- ComputeRSSPhase(lc,omega_est,tem,phis=phis)
plot(phis,rss_phi)
abline(v=coeffs[4])



##
## PART 2: refit model, subtracting off dust first
##         here we test if model works when use.dust=FALSE
ebv <- coeffs[2]
## remove dust
bands <- names(tem$dust)
for(ii in 1:length(bands)){
    lc$mag[lc$band==bands[ii]] <- lc$mag[lc$band==bands[ii]] - tem$dust[ii]*ebv
}


rss <- FitTemplate(lc,omegas,tem,use.dust=FALSE)
omega_est <- omegas[which.min(rss)]
p_est <- 1/omega_est
coeffs <- ComputeCoeffs(lc,omega_est,tem,use.dust=FALSE)
names(coeffs) <- c("mu","ebv","amp","phase")

## view rss
plot(1/omegas,rss,xlab="period",ylab="rss")
abline(v=p_est)

## plot folded light curve with best fit
plotLC(lc,p_est,coeffs,tem)



## check that phase determine by NewtonUpdate is close to best phase
## by performing a grid search
phis <- (1:100)/100
rss_phi <- ComputeRSSPhase(lc,omega_est,tem,phis=phis,use.dust=FALSE)
plot(phis,rss_phi)
abline(v=coeffs[4])



##
## PART 3: make some observations in lc very bad
##         
ix <- sample(1:nrow(lc),floor(nrow(lc)*.6))
lc[ix,3] <- lc[ix,3] + rnorm(length(ix),mean=0,sd=3)
lc[ix,4] <- sqrt(lc[ix,4]^2 + 3^2)

## plot raw light curve
colpch <- 1:5
names(colpch) <- unique(lc$band)
plot(lc$time,lc$mag,col=colpch[lc$band],pch=colpch[lc$band],
     ylim=rev(range(lc$mag)),
     xlab="time",ylab="magnitude")
segments(lc$time,lc$mag+lc$error,lc$time,lc$mag-lc$error)


###### a) use.errors=TRUE
rss <- FitTemplate(lc,omegas,tem)
omega_est <- omegas[which.min(rss)]
p_est <- 1/omega_est
coeffs <- ComputeCoeffs(lc,omega_est,tem)
names(coeffs) <- c("mu","ebv","amp","phase")

## view rss
plot(1/omegas,rss,xlab="period",ylab="rss")
abline(v=p_est)

## plot folded light curve with best fit
plotLC(lc,p_est,coeffs,tem)

## check that phase determine by NewtonUpdate is close to best phase
## by performing a grid search
phis <- (1:100)/100
rss_phi <- ComputeRSSPhase(lc,omega_est,tem,phis=phis)
plot(phis,rss_phi)
abline(v=coeffs[4])


###### b) use.errors=FALSE, should get poor parameter estimates
rss <- FitTemplate(lc,omegas,tem,use.errors=FALSE)
omega_est <- omegas[which.min(rss)]
p_est <- 1/omega_est
coeffs <- ComputeCoeffs(lc,omega_est,tem,use.errors=FALSE)
names(coeffs) <- c("mu","ebv","amp","phase")

## view rss
plot(1/omegas,rss,xlab="period",ylab="rss")
abline(v=p_est)

## plot folded light curve with best fit
plotLC(lc,p_est,coeffs,tem)

## check that phase determine by NewtonUpdate is close to best phase
## by performing a grid search
phis <- (1:100)/100
rss_phi <- ComputeRSSPhase(lc,omega_est,tem,phis=phis,use.errors=FALSE)
plot(phis,rss_phi)
abline(v=coeffs[4])

