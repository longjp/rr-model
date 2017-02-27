## this demo shows how sinusoidal
## period estimation algorithms work
## by James Long

rm(list=ls())

construct_design <- function(w,K,t){
    predesign <- w*outer(t,1:K)
    return(cbind(1,cos(predesign),sin(predesign)))
}

get_freqs <-function(period_min,period_max,freq_del = 0.1/4000){
    freq_max <- 1/period_min
    freq_min <- 1/period_max
    return(2 * pi * seq(freq_min, freq_max, freq_del))
}

compute_params <- function(w,K,mag,weights,X){
    B <- t(X) %*% (X * weights)
    d <- t(X) %*% (mag * weights)
    return(solve(B,d))
}   

compute_rss <- function(w,K,lc){
    X <- construct_design(w,K,lc[,1])
    beta <- compute_params(w,K,lc[,2],lc[,3]^{-2},X)
    r <- (lc[,2] - X%*%beta)
    return(sum(lc[,3] * (r^2)))
}


## lc should be matrix or dataframe with 3 columns
## (times,magnitudes,errors)
lc <- read.table("AllLCs/LC_4099.dat")
K <- 2 ## try other K as well
period_min <- 0.2 ## since this is a cepheid, should have period in [1,100]
period_max <- 1
omegas <- get_freqs(period_min,period_max,.1/diff(range(lc[,1])))

bands <- unique(lc[,2])
bands <- sort(as.character(bands))
lcs <- list()
for(ii in 1:length(bands)){
    lcs[[ii]] <- lc[lc[,2]==bands[ii],c(1,3,4)]
}
names(lcs) <- bands

rss <- vapply(lcs,function(x){vapply(omegas,compute_rss,c(0),K,x)},rep(0,length(omegas)))
rss <- rowSums(rss)
plot(omegas,rss)
abline(v=omegas[which.min(rss)])

## the omega estimate is the omega with the lowest rss 
p_est <- (2*pi)/omegas[which.min(rss)]



cex.lab <- 2
cex.axis <- 2
err.scale <- 2


xlim <- range(lc[,1])
ylim <- range(lc[,3])
pdf("../figs/rrlyrae_nomodel_fit.pdf",height=6,width=12)
par(mar=c(5,5,1,1))
plot(0,0,xlim=xlim,ylim=rev(ylim),xaxs='i',xlab="Time",ylab="Magnitude",
     cex.lab=cex.lab,cex.axis=cex.axis)
for(ii in 1:length(bands)){
    segments(lcs[[ii]][,1],lcs[[ii]][,2]-err.scale*lcs[[ii]][,3],
             lcs[[ii]][,1],lcs[[ii]][,2]+err.scale*lcs[[ii]][,3],col='grey')
    points(lcs[[ii]][,1],lcs[[ii]][,2],ylim=rev(range(lc[,3])),col=ii,pch=ii)
}
nband <- length(bands)
legend("bottomleft",paste0(bands," Band"),col=1:nband,pch=1:nband,cex=1.3)
dev.off()


######### sine with 2 harmonics function
## draw model on plot
xlim <- c(0,1)
ylim <- range(lc[,3])
pdf("../figs/rrlyrae_model_fit.pdf",height=6,width=12)
par(mar=c(5,5,1,1))
plot(0,0,xlim=xlim,ylim=rev(ylim),xaxs='i',
     xlab=paste0("Phase (period=",round(p_est,2),")"),ylab="Magnitude",cex.lab=cex.lab,cex.axis=cex.axis)
for(ii in 1:length(bands)){
    points((lcs[[ii]][,1]%%p_est)/p_est,lcs[[ii]][,2],col=ii,pch=ii)
    segments((lcs[[ii]][,1]%%p_est)/p_est,lcs[[ii]][,2] + err.scale*lcs[[ii]][,3],
    (lcs[[ii]][,1]%%p_est)/p_est,lcs[[ii]][,2] - err.scale*lcs[[ii]][,3],
    col='grey')
    X <- construct_design(2*pi/p_est,K,lcs[[ii]][,1])
    beta <- compute_params(2*pi/p_est,K,lcs[[ii]][,2],lcs[[ii]][,3]^2,X)
    t <- (lcs[[ii]][,1] %% p_est)/p_est
    m <- X%*%beta
    m <- m[order(t)]
    t <- t[order(t)]
    points(t,m,col=ii,lty=ii,pch=19,type='l')
}
legend("bottomleft",paste0(bands," Band"),col=1:nband,pch=1:nband,cex=1.3)
dev.off()






######### pure sine function
## draw model on plot
K <- 1
xlim <- c(0,1)
ylim <- range(lc[,3])
pdf("../figs/rrlyrae_model_fit_sine.pdf",height=6,width=12)
par(mar=c(5,5,1,1))
plot(0,0,xlim=xlim,ylim=rev(ylim),xaxs='i',
     xlab=paste0("Phase (period=",round(p_est,2),")"),ylab="Magnitude",cex.lab=cex.lab,cex.axis=cex.axis)
for(ii in 1:length(bands)){
    points((lcs[[ii]][,1]%%p_est)/p_est,lcs[[ii]][,2],col=ii,pch=ii)
    segments((lcs[[ii]][,1]%%p_est)/p_est,lcs[[ii]][,2] + err.scale*lcs[[ii]][,3],
    (lcs[[ii]][,1]%%p_est)/p_est,lcs[[ii]][,2] - err.scale*lcs[[ii]][,3],
    col='grey')
    X <- construct_design(2*pi/p_est,K,lcs[[ii]][,1])
    beta <- compute_params(2*pi/p_est,K,lcs[[ii]][,2],lcs[[ii]][,3]^2,X)
    t <- (lcs[[ii]][,1] %% p_est)/p_est
    m <- X%*%beta
    m <- m[order(t)]
    t <- t[order(t)]
    points(t,m,col=ii,lty=ii,pch=19,type='l')
}
legend("bottomleft",paste0(bands," Band"),col=1:nband,pch=1:nband,cex=1.3)
dev.off()



