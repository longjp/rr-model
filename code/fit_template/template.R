rm(list=ls())
source("fit_template.R")
load("template.RData")

## read in data and plot
fname <- "LC_4099.dat"
lc <- read.table(fname)
names(lc) <- c("time","band","mag","error")
lc$band <- as.character(lc$band)
lc$time <- lc$time - min(lc$time)
## plot raw light curve
colpch <- 1:5
names(colpch) <- unique(lc$band)
plot(lc$time,lc$mag,col=colpch[lc$band],pch=colpch[lc$band],
     ylim=rev(range(lc$mag)),
     xlab="time",ylab="magnitude")
segments(lc$time,lc$mag+lc$error,lc$time,lc$mag-lc$error)


## fit template model and obtain coefficients
omegas <- seq(from=1.0,to=5.0,by=0.1/4000.0)
rss <- FitTemplate(lc,omegas,tem)
omega_est <- omegas[which.min(rss)]
p_est <- 1/omega_est
coeffs <- ComputeCoeffs(lc,omega_est,tem)
names(coeffs) <- c("mu","ebv","amp","phase")

## view rss
plot(1/omegas,rss,xlab="period",ylab="rss")
abline(v=p_est)

## plot folded light curve with best fit
colpch <- 1:5
names(colpch) <- unique(lc$band)
plot((lc$time %% p_est)/p_est,lc$mag,
     col=colpch[lc$band],pch=colpch[lc$band],
     ylim=rev(range(lc$mag)),
     xlab="time",ylab="magnitude",
     xlim=c(0,1),xaxs='i')
segments((lc$time %% p_est)/p_est,
         lc$mag+lc$error,
         (lc$time %% p_est)/p_est,
         lc$mag-lc$error)
ti <- (1:100)/100
m <- PredictAllBand(ti,1,coeffs,tem)
for(ii in 1:length(tem$betas)){
    points(ti,m[,ii],type='l',col=colpch[names(tem$betas)[ii]])
}
