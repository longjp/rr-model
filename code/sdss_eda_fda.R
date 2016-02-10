rm(list=ls())

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


### extract location max, location min, amp, beta0 for each lc,band
params <- array(0,dim=c(length(tmss),5,4),dimnames=list(NULL,bands,c("max","min","amp","beta")))
out <- list()
for(ii in 1:length(tmss)){
    out[[ii]] <- list()
    for(jj in 1:length(bands)){
        lc <- tmss[[ii]][[bands[jj]]]
        lc[,1] <- (lc[,1] %% rrlyrae[ii,3]) / rrlyrae[ii,3]
        temp <- supsmu(lc[,1],lc[,2],periodic=TRUE)
        out[[ii]][[jj]] <- temp
        params[ii,jj,] <- c(temp$x[which.max(temp$y)],temp$x[which.min(temp$y)],max(temp$y)-min(temp$y),mean(temp$y))
    }
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


## find cc
phis <- (params[,,"min"] - params[,,"max"]) %% 1
phis <- as.vector(phis)
phis <- cbind(cos(2*pi*phis),sin(2*pi*phis))
phis <- colSums(phis)
phis_av <- phis / sqrt(sum(phis^2))
cc <- atan2(phis_av[2],phis_av[1]) / (2*pi)


######## find phij
phis <- rep(0,length(bands))
names(phis) <- bands
for(jj in 1:length(bands)){
    phis_temp <- (params[,bands[jj],"min"] - params[,1,"min"]) %% 1
    phis_temp <- cbind(cos(2*pi*phis_temp),sin(2*pi*phis_temp))
    phis_temp <- colSums(phis_temp)
    phis_av <- phis_temp / sqrt(sum(phis_temp^2))
    phis[jj] <- atan2(phis_av[2],phis_av[1]) / (2*pi)
}


save(betas,amps,phis,cc,file="sdss_eda.RData")
