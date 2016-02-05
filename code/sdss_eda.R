rm(list=ls())

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
for(ii in 1:length(tmss)){
    for(jj in 1:length(bands)){
        lc <- tmss[[ii]][[bands[jj]]]
        lc[,1] <- (lc[,1] %% rrlyrae[ii,3]) / rrlyrae[ii,3]
        out <- supsmu(lc[,1],lc[,2],periodic=TRUE)
        params[ii,jj,] <- c(out$x[which.max(out$y)],out$x[which.min(out$y)],max(out$y)-min(out$y),mean(out$y))
    }
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
    betas[jj] <- median(params[,jj,"beta"] - params[,1,"beta"])
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
    phis_temp <- (params[,jj,"min"] - params[,1,"min"]) %% 1
    phis_temp <- cbind(cos(2*pi*phis_temp),sin(2*pi*phis_temp))
    phis_temp <- colSums(phis_temp)
    phis_av <- phis_temp / sqrt(sum(phis_temp^2))
    phis[jj] <- atan2(phis_av[2],phis_av[1]) / (2*pi)
}


save(betas,amps,phis,cc,file="sdss_eda.RData")
