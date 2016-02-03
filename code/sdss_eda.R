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
    print(bands[ii])
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
        plot(lc[,1],lc[,2])
        out <- supsmu(lc[,1],lc[,2],periodic=TRUE)
        params[ii,jj,] <- c(out$x[which.max(out$y)],out$x[which.min(out$y)],max(out$y)-min(out$y),mean(out$y))
    }
}



pairs(params[,,"amp"])

pairs(params[,,"beta"])



for(jj in 2:5){
    print(jj)
    print(lm(params[,1,"beta"] ~ params[,jj,"beta"]))
}




for(jj in 2:5){
    print(jj)
    print(lm(params[,1,"amp"] ~ params[,jj,"amp"]))
}



pairs(params[,,2]-params[,,1])



#### issues:
## 1. what band to define as having ampj=1 and phij = 0
## 2. deal with magnitudes being inversely prop to brightness (should we flip RR model)
