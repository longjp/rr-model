rm(list=ls())

source('func_sine.R')
source('func.R')

library(RColorBrewer)
library(scales)
library(fda)

rrlyrae <- read.table("apj326724t2_mrt.txt",skip=42)
rrlyrae <- rrlyrae[rrlyrae[,2] == "ab",]

folder <- "rrlyrae"
fnames <- list.files(folder)



## order fnames and rrlyrae
rrlyrae[,1] <- paste("LC_",rrlyrae[,1],".dat",sep="")
fnames <- fnames[fnames %in% rrlyrae[,1]]
rrlyrae <- rrlyrae[rrlyrae[,1] %in% fnames,]
rrlyrae <- rrlyrae[order(rrlyrae[,1]),]
fnames <- fnames[order(fnames)]


## load light curves and put in tmss format
lcs <- list()
for(ii in 1:length(fnames)){
    lcs[[ii]] <- read.table(paste(folder,fnames[ii],sep="/"))
}
for(ii in 1:length(lcs)) names(lcs[[ii]]) <- c("time","b","mag","sigma")


LCtoTMS <- function(lc){
    lc[,1] <- lc[,1] - min(lc[,1])
    levs <- levels(lc$b)
    tms <- list()
    for(ii in 1:length(levs)){
        tms[[ii]] <- lc[lc$b == levs[ii],c("time","mag","sigma")]
        names(tms[[ii]]) <- c("time","mag","sigma")
    }
    names(tms) <- levs
    return(tms)
}
tmss <- lapply(lcs,LCtoTMS)



## get light curves with at least 50 observations / band in each band
bands <- names(tmss[[1]])
nobs <- matrix(0,nrow=length(tmss),ncol=5)
for(ii in 1:length(bands)){
    nobs[,ii] <- vapply(tmss,function(x){nrow(x[[bands[ii]]])},c(0))
}
to_use <- rowSums(nobs > 45) == 5
tmss <- tmss[to_use]
rrlyrae <- rrlyrae[to_use,]

N <- 100
t <- (1:N)/N
### extract location max, location min, amp, beta0 for each lc,band
params <- array(0,dim=c(length(tmss),5,4),dimnames=list(NULL,bands,c("max","min","amp","beta")))
out <- list()
lc_grid <- array(0,c(length(tmss),N,5))
for(ii in 1:length(tmss)){
    out[[ii]] <- list()
    for(jj in 1:length(bands)){
        lc <- tmss[[ii]][[bands[jj]]]
        lc[,1] <- (lc[,1] %% rrlyrae[ii,3]) / rrlyrae[ii,3]
        temp <- supsmu(lc[,1],lc[,2],periodic=TRUE)
        ymean <- mean(c(temp$y[1],temp$y[length(temp$y)]))
        temp$y <- c(ymean,temp$y,ymean)
        temp$x <- c(0,temp$x,1)
        out[[ii]][[jj]] <- temp
        lc_grid[ii,,jj] <- approx(temp$x,temp$y,xout=t,rule=2)$y
        params[ii,jj,] <- c(temp$x[which.max(temp$y)],temp$x[which.min(temp$y)],max(temp$y)-min(temp$y),mean(temp$y))
    }
}



alphas <- apply(lc_grid,1,mean)
lc_band_means <- apply(lc_grid,c(1,3),mean)
betas <- colMeans(lc_band_means - alphas)


colnames(lc_band_means) <- bands
pairs(lc_band_means)





cols <- brewer.pal(10,name="RdBu")
decLocations <- quantile(rrlyrae[,3], probs = seq(0.1,0.9,by=0.1),type=4)
dec <- findInterval(rrlyrae[,3],c(-Inf,decLocations, Inf))


lc_ix <- as.factor(rep(1:nrow(lc_band_means),ncol(lc_band_means)))
bands_ix <- as.factor(rep(1:ncol(lc_band_means),each=nrow(lc_band_means)))
lc_means <- as.vector(lc_band_means)

lm.fit <- lm(lc_means ~ lc_ix + bands_ix)
resid <- matrix(lm.fit$residuals,ncol=length(bands))
colnames(resid) <- bands
pairs(resid,col=cols[dec])


resid_svd <- svd(resid)
pred <- resid_svd$d[1]*resid_svd$u[,1,drop=FALSE]%*%matrix(resid_svd$v[,1],ncol=5)
dev.new()
pairs(resid - pred,col=cols[dec])


pairs(cbind(resid - pred,log(rrlyrae[,3])),col=cols[dec])

resid2 <- resid - pred

## regress residuals on period and get residuals of residuals
for(ii in 1:ncol(resid2)){
    resid2[,ii] <- lm(resid2[,ii] ~ rrlyrae[,3])$residuals
}

dev.new()
pairs(cbind(resid2,rrlyrae[,3]),col=cols[dec])






#### CORRELATION IN RESIDUALS FOR G AND U APPEARS RELATED TO PERIOD





for(ii in 1:length(tmss)){
    lc_grid[ii,,] <- lc_grid[ii,,] - alphas[ii]
    for(jj in 1:5){
        lc_grid[ii,,jj] - betas[jj]
    }
}


## phase and amplitude align light curves in lc_grid
for(ii in 1:length(tmss)){
    temp <- lc_grid[ii,,1]
    ##temp <- temp - mean(temp)
    ##temp <- temp / (max(temp) - min(temp))
    ix <- which.min(temp)
    if(ix > 1.5){
        for(jj in 1:5){
            temp <- lc_grid[ii,,jj]
            lc_grid[ii,,jj] <- (c(temp[ix:N],temp[1:(ix-1)]) - mean(temp))/(max(temp)-min(temp))
        }
    } else {
        for(jj in 1:5){
            temp <- lc_grid[ii,,jj]
            lc_grid[ii,,jj] <- (temp - mean(temp))/(max(temp)-min(temp))
        }
    }
}

cols <- brewer.pal(10,name="RdBu")
decLocations <- quantile(rrlyrae[,3], probs = seq(0.1,0.9,by=0.1),type=4)
dec <- findInterval(rrlyrae[,3],c(-Inf,decLocations, Inf))



del <- phase(lc_grid[,,1])
for(JJ in 1:5){
    x_new <- list()
    n <- nrow(lc_grid[,,JJ])
    x_new[[JJ]] <- t(vapply(1:n,function(ii){phase_shift(lc_grid[ii,,JJ],del[ii])},rep(0,N)))
    dev.new()
    par(mfcol=c(2,1))
    ylim <- range(lc_grid[,,JJ])
    plot(0,0,ylim=ylim,xlim=c(0,1),col=0,xlab="phase",ylab="mag",xaxs="i")
    for(ii in 1:length(tmss)){
        points(t,lc_grid[ii,,JJ],type='l',col=cols[dec[ii]])
    }
    ylim <- range(x_new[[JJ]])
    plot(0,0,col=0,ylim=ylim,xlim=c(0,1),xlab="",ylab="",xaxs="i")
    for(ii in 1:nrow(x_new[[JJ]])){
        points(t,x_new[[JJ]][ii,],type='l',col=cols[dec[ii]])
    }
}











JJ <- 2
n <- nrow(lc_grid[,,JJ])
x_new[[JJ]] <- t(vapply(1:n,function(ii){phase_shift(lc_grid[ii,,JJ],del[ii])},rep(0,N)))
dev.new()
par(mfcol=c(2,1))
ylim <- range(lc_grid[,,JJ])
plot(0,0,ylim=ylim,xlim=c(0,1),col=0,xlab="phase",ylab="mag",xaxs="i")
for(ii in 1:length(tmss)){
    points(t,lc_grid[ii,,JJ],type='l',col="#00000030")
}
ylim <- range(x_new[[JJ]])
plot(0,0,col=0,ylim=ylim,xlim=c(0,1),xlab="",ylab="",xaxs="i")
for(ii in 1:nrow(x_new[[JJ]])){
    points(t,x_new[[JJ]][ii,],type='l',col="#00000030")
}











## construct lightcurve model from out
for(ii in 1:length(out)){
    phase <- out[[ii]][[1]]$x[which.min(out[[ii]][[1]]$y)]
    for(jj in 1:5){
        out[[ii]][[jj]]$x <- (out[[ii]][[jj]]$x - phase) %% 1
        out[[ii]][[jj]]$y <- (out[[ii]][[jj]]$y - params[ii,jj,4]) / params[ii,jj,3]
    }
}

ylim <- c(-0.75,.75)
plot(0,0,xlim=c(0,1),ylim=ylim,col=0,xlab="phase",ylab="mag")
for(ii in 1:length(out)){
    for(jj in 1:5){
        ords <- order(out[[ii]][[jj]]$x)
        points(out[[ii]][[jj]]$x[ords],out[[ii]][[jj]]$y[ords],col=alpha(jj, 0.1),type='l')
    }
}


cols <- brewer.pal(10,name="RdBu")
decLocations <- quantile(rrlyrae[,3], probs = seq(0.1,0.9,by=0.1),type=4)
dec <- findInterval(rrlyrae[,3],c(-Inf,decLocations, Inf))


dev.new()
JJ <- 1
ylim <- c(-0.75,.75)
plot(0,0,xlim=c(0,1),ylim=ylim,col=0,xlab="phase",ylab="mag")
for(ii in 1:length(out)){
    ords <- order(out[[ii]][[JJ]]$x)
    points(out[[ii]][[JJ]]$x[ords],out[[ii]][[JJ]]$y[ords],type='l',col=cols[dec[ii]])
}



dev.new()
JJ <- 2
ylim <- c(-0.75,.75)
plot(0,0,xlim=c(0,1),ylim=ylim,col=0,xlab="phase",ylab="mag")
for(ii in 1:length(out)){
    ords <- order(out[[ii]][[JJ]]$x)
    points(out[[ii]][[JJ]]$x[ords],out[[ii]][[JJ]]$y[ords],col=cols[dec[ii]],type='l')
}


dev.new()
JJ <- 3
ylim <- c(-0.75,.75)
plot(0,0,xlim=c(0,1),ylim=ylim,col=0,xlab="phase",ylab="mag")
for(ii in 1:length(out)){
        ords <- order(out[[ii]][[JJ]]$x)
        points(out[[ii]][[JJ]]$x[ords],out[[ii]][[JJ]]$y[ords],col=alpha(JJ, 0.1),type='l')
}


dev.new()
JJ <- 4
ylim <- c(-0.75,.75)
plot(0,0,xlim=c(0,1),ylim=ylim,col=0,xlab="phase",ylab="mag")
for(ii in 1:length(out)){
        ords <- order(out[[ii]][[JJ]]$x)
        points(out[[ii]][[JJ]]$x[ords],out[[ii]][[JJ]]$y[ords],col=alpha(JJ, 0.1),type='l')
}


dev.new()
JJ <- 5
ylim <- c(-0.75,.75)
plot(0,0,xlim=c(0,1),ylim=ylim,col=0,xlab="phase",ylab="mag")
for(ii in 1:length(out)){
        ords <- order(out[[ii]][[JJ]]$x)
        points(out[[ii]][[JJ]]$x[ords],out[[ii]][[JJ]]$y[ords],col=alpha(JJ, 0.1),type='l')
}





### get ratio of amplitudes, g band is defined as 0
amps <- rep(0,length(bands))
names(amps) <- bands
for(jj in 1:length(bands)){
    amps[jj] <- lm(params[,bands[jj],"amp"] ~ params[,bands[1],"amp"]-1)$coefficients
}

### get median difference in mean mags across bands
betas <- rep(0,length(bands))
names(betas) <- bands
for(jj in 1:length(bands)){
    betas[jj] <- median(params[,bands[jj],"beta"] - params[,1,"beta"])
}




save(betas,amps,phis,cc,file="sdss_eda.RData")
