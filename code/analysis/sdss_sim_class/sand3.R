### created for cdi grant proposal
### downsample RRL multiple times, compute rss curve
### are rss local minimia always in same region of frequency space
rm(list=ls())

figfold <- function(figname){
    paste0("sand3/",figname)
}

## load necessary libraries
load("../../fit_template/template_sdss.RData")
source("../../fit_template/fit_template.R")
source("../../common/funcs.R")
source("../funcs.R")
source('sine_funcs.R')

## data source
load("../../data/clean/sdss_sim_class.RData")

## parameters for simulation
source("../params.R")

## create dust corrected tms
ebv <- extcr / tem$dust['r']
tmsc <- mapply(DustCorrect,tms,ebv,MoreArgs=list(tem=tem),SIMPLIFY=FALSE)
tmsc_FULL <- mapply(DustCorrect,tms_FULL,ebv,MoreArgs=list(tem=tem),SIMPLIFY=FALSE)


ix <- which.max(vapply(tmsc_FULL[1:sum(cl=="rr")],function(x){sum(vapply(x,nrow,c(0)))},c(0)))
ix
lapply(tmsc_FULL[[ix]],nrow)

tm <- tmsc_FULL[[ix]][['g']]
tm <- tm[tm[,3] < 0.02,]
p <- periods[ix]
tm



png(figfold("data.png"),width=800,height=400)
plot(tm[,1],tm[,2],xlab="t",ylab="y")
dev.off()


omegas <- get_freqs(0.4,0.9)
t <- tm[,1]
y <- tm[,2]
rss <- vapply(omegas,function(w){compute_rss(w,1,t,y)},c(0))
omega_true <- omegas[which.min(rss)]

png(figfold("rss.png"),width=800,height=400)
plot(omegas,rss,xlim=range(omegas),ylim=range(rss),xlab="Frequency",ylab="RSS",col="#00000030")
abline(v=omega_true,col='red',lwd=2)
dev.off()

omega_true
2*pi/omega_true


png(figfold("data_folded.png"),width=800,height=400)
plot(tm[,1] %% p,tm[,2],xlab="Phase = t modulo period",ylab="y")
dev.off()






N <- 20
rss <- matrix(0,nrow=N,ncol=length(omegas))
omega <- rep(0,N)
for(ii in 1:N){
    print(ii)
    ix <- sample(1:length(y),15)
    rss[ii,] <- vapply(omegas,function(w){compute_rss(w,1,t[ix],y[ix])},c(0))
    omega[ii] <- omegas[which.min(rss[ii,])]
}

rss <- rss - apply(rss,1,max)
worst_true <- max(rss[,omega_true])


png(figfold("rss_down.png"),width=800,height=400)
plot(0,0,xlim=range(omegas),ylim=range(rss),col=0,xlab="Frequency",ylab="RSS")
abline(v=omega,col='yellow',lwd=2)
for(ii in 1:N){
    points(omegas,rss[ii,],col='#00000010')
}
##rmed <- apply(rss,2,function(x){quantile(x,0.1)})
rmed <- apply(rss,2,median)
points(omegas,rmed,type='l',lwd=2,col='red')
##points(omega_true,worst_true,lwd=2,pch=3,col='yellow')
dev.off()



ix <- which.min(omega)
ix

png(figfold("rss_one_down.png"),width=800,height=400)
plot(omegas,rss[ix,],xlim=range(omegas),xlab="Frequency",ylab="RSS",col="#00000050")
abline(v=omega[ix],col='yellow',lwd=2)
abline(v=omega_true,col='red',lwd=2)
dev.off()

## the max rss varies a lot by run

## abline(v=2*pi/p,col='blue',lwd=2)
## abline(v=omega,col='red',lwd=2)

