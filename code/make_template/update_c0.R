rm(list=ls())
unlink("*.RData")
source('../common/funcs.R')
source('plot_options.R')
source('../fit_template/fit_template.R')
load("../data/clean/sdss_rrab.RData")
library(ellipse)

plot_foldername <- "figs"

## remove photometric measurements with uncertainty greater than scut
scut <- .2
for(ii in 1:length(tms)){
    for(jj in 1:length(tms[[ii]])){
        temp <- tms[[ii]][[jj]]
        tms[[ii]][[jj]] <- temp[temp[,3] < scut,]
    }
}    

## get light curves with at least 50 observations / band in each band
bands <- names(tms[[1]])
bands <- bands[order(bands)]
nobs <- matrix(0,nrow=length(tms),ncol=length(bands))
for(ii in 1:length(bands)){
    nobs[,ii] <- vapply(tms,function(x){nrow(x[[bands[ii]]])},c(0))
}
to_use <- rowSums(nobs > 45) == 5
tms <- tms[to_use]
periods <- periods[to_use]
Nlc <- length(tms)

## interpolate light curves on a grid
N <- 100
t <- (1:N)/N
lc_grid <- array(0,c(Nlc,N,5),dimnames=list(NULL,NULL,bands))
for(ii in 1:Nlc){
    for(jj in 1:length(bands)){
        lc <- tms[[ii]][[bands[jj]]]
        lc[,1] <- (lc[,1] %% periods[ii]) / periods[ii]
        ords <- order(lc[,1])
        x <- lc[ords,1]
        y <- lc[ords,2]
        ymean <- mean(c(y[1],y[length(y)]))
        y <- c(ymean,y,ymean)
        x <- c(0,x,1)
        y_approx <- approx(x,y,xout=t,rule=2)$y
        lc_grid[ii,,jj] <- y_approx
    }
}



## ## compute mean mag in each band / lc
m <- apply(lc_grid,c(1,3),mean)
pdf("figs/band_means.pdf")
pairs(m)
dev.off()

## load extinction parameters
dat <- read.table("extc.dat")
dat <- dat[dat$V1=="SDSS",c(2,3)]
dust <- dat[,2]
names(dust) <- dat[,1]
dust <- dust[order(names(dust))]

## load and compute beta (M_x) parameters using median of periods
rrmag <- read.table("rrmag.dat",header=TRUE)
rrmag <- rrmag[rrmag$Sys=="SDSS",]
rrmag <- rrmag[order(rrmag$bnd),]

## absolute mag period dependence
betasM <- matrix(0,nrow=length(periods),ncol=5)
for(ii in 1:length(periods)){
    betasM[ii,] <- rrmag$c0 + rrmag$p1*(log10(periods[ii]) + 0.2) + rrmag$p2*(log10(periods[ii]) + 0.2)^2
}
colnames(betasM) <- rrmag$bnd


## m <- matrix(0,nrow=length(tms),ncol=length(bands))
## for(ii in 1:length(bands)){
##     m[,ii] <- vapply(tms,function(x){mean(x[[bands[ii]]][,2])},c(0))
## }



## estimate ebv and mu for each lc
ests <- data.frame(ID=names(tms),mu_est=0,extcr_est=0)

## compare dust to sesar dust
m_shift <- m - betasM
rss <- rep(0,nrow(m))
for(ii in 1:nrow(m)){
    lm_fit <- lm(m_shift[ii,]~dust)
    ests[ii,2:3] <- lm_fit$coefficients
    rss[ii] <- sum(lm_fit$residuals^2)/(length(dust)-2)
}
ests[,3] <- ests[,3]*dust['r']



kpc_to_mu <- function(kpc) 5*(log10(kpc*1000)-1)
mu_to_kpc <- function(mu) (10^(mu/5 + 1)) / 1000

tab <- read.table("../data/raw/apj326724t3_mrt.txt",skip=30)
tab <- tab[,c(1,4,5)]
tab[,3] <- kpc_to_mu(tab[,3])
names(tab) <- c("ID","extcr","mu")
tab$ID <- paste0("LC_",tab$ID,".dat")

comp <- merge(tab,ests)
plot(comp$mu,comp$mu_est)
plot(comp$extcr,comp$extcr_est)


y <- comp$mu_est - comp$mu
x <- comp$extcr_est-comp$extcr
par(mar=c(5,5,1,1))
plot(x,y,xlab="Extinction r Residual",ylab="mu Residual")
lm_fit <- lm(y ~ x)
abline(lm_fit$coeff)
abline(v=0)
abline(h=0)

X <- cbind(dust/dust['r'],1)
covx <- solve(t(X)%*%X)
covx <- (var(x)/covx[1,1])*covx
##covx <- solve(t(X)%*%X) * mean(rss) analytic method

mean(x)

points(ellipse(covx,centre=c(mean(x),mean(y))),col='red',type='l',pch=2,lwd=2)
points(mean(x),mean(y),pch=19,col='blue',cex=1.5)



## update c0
cc <- mean(x) / dust['r']
c0 <- rrmag$c0 + cc*dust


save(c0,file="update_c0.RData")
