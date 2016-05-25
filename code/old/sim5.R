### compare newton method to grid search
### for newton, condition on amplitude, beta0, update phase
### optimize only across phase

rm(list=ls())
library(scales)
library(parallel)
source("func_saw_multi.R")
source("func_multi.R")

omegas <- seq(1,5,length.out=40000)

## construct data with panstarrs cadence
# parameters
phi <- runif(1)
amp <- 0.3
beta0 <- rnorm(1)
omega <- sample(omegas,1)
ix <- 9


## ## panstarrs cadence
f <- list.files("PS1_sample_LCs",full.names=TRUE)
tab <- read.table(f[ix],header=TRUE)
t <- tab[,3]
e.sd <- tab[,1]
band <- tab[,2]

## simulate times, e.d, band
## n <- 100
## t <- cumsum(rexp(n))
## e.sd <- rep(.3,n)
## band <- sample(c("red","blue"),n,replace=TRUE)


## construct amplitude and phase offsets for bands
B <- length(unique(band))
amp_mag <- (0:(B-1))/10 + 1
names(amp_mag) <- unique(band)
ampj <- amp_mag[band]
phi_mag <- (0:(B-1))/10
names(phi_mag) <- unique(band)
phij <- phi_mag[band]
## 
m <- amp*ampj*Saw(omega*t+phij+phi) + beta0 + rnorm(length(e.sd),mean=0,sd=e.sd)
plot(t,m)


plot(t %% (1/omega),m,col=band)




#### grid search to estimate omega and phase, closed for solution for beta0 and amp
N <- 500
phis_grid <- (0:(N-1))/N
mc.cores <- 7
tm <- proc.time()
## rss <- rep(0,length(omegas))
## for(ii in 1:length(rss)){
##     print(ii)
##     ## m.resid <- vapply(phis_grid,function(x){ComputeResiduals(m,t,x,omegas[ii])},rep(0,length(t))) ## faster not using lm
##     ## rss[ii] <- min(colSums(m.resid^2))
##     m.resid <- mclapply(phis_grid,function(x){ComputeResiduals(m,t,x,omegas[ii])},mc.cores=mc.cores) ## faster not using lm
##     rss[ii] <- min(vapply(m.resid,function(x){sum(x^2)},c(0)))
## }
compute_rss_parallel <- function(ii){
    print(ii)
    m.resid <- vapply(phis_grid,function(x){ComputeResiduals(m,t,ampj,phij,x,omegas[ii])},rep(0,length(t))) ## faster not using lm
    return(min(colSums(m.resid^2)))
}
rss <- unlist(mclapply(1:length(omegas),function(ii){compute_rss_parallel(ii)},mc.cores=mc.cores))
proc.time() - tm
plot(omegas,rss)
abline(v=omega,col="black")
abline(v=omegas[which.min(rss)],col="blue")

## use newton method
tm <- proc.time()
set.seed(1234)
rss_newton <- SawRss(m,t,ampj,phij,omegas,NN=1)
proc.time() - tm
points(omegas,rss_newton,col='red')
abline(v=omegas[which.min(rss_newton)],col="red",lwd=2)


save(m,t,omega,file="sim4_newton_fail_2.RData")



