## appears to be period evolution of lc shape for rrab
rm(list=ls())
library(splines)
library(KernSmooth)



######## 1. OBTAIN AVERAGE RR LYRAE AMPLITUDE ACROSS ALL BANDS
stars <- read.table("apj326724t2_mrt.txt",skip=42)
folder <- "rrlyrae"
fnames <- list.files(folder)


## put fnames and stars in order
stars$V1 <- paste0("LC_",as.character(stars$V1),".dat")
stars <- stars[stars$V1 %in% fnames,]
stars <- stars[stars[,2] == "ab",]

periods <- stars[,3]
fnames <- stars[,1]

## read in files and convert to standard format
lcs <- list()
for(ii in 1:length(fnames)){
    lcs[[ii]] <- read.table(paste(folder,fnames[ii],sep="/"))
}
for(ii in 1:length(lcs)) names(lcs[[ii]]) <- c("time","b","mag","sigma")

ConvertSDSS <- function(lc){
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
lcs <- lapply(lcs,ConvertSDSS)

for(ii in 1:length(lcs[[1]])){
    print(names(lcs[[1]])[[ii]])
    print(summary(vapply(lcs,function(x){nrow(x[[ii]])},c(0))))
}


## get rid of everything with very few observations
n_points <- vapply(lcs,function(x){min(vapply(x,function(y){nrow(y)},c(0)))},c(0))
to_use <- n_points > 45
periods <- periods[to_use]
lcs <- lcs[to_use]





## sample at N time points
N <- 200
sm <- array(0,c(length(lcs),length(lcs[[1]]),N))
amps <- matrix(0,length(lcs),length(lcs[[1]]))
for(ii in 1:length(lcs)){
    for(jj in 1:length(lcs[[ii]])){
        temp <- lcs[[ii]][[jj]]
        temp[,1] <- (temp[,1] %% periods[ii]) / periods[ii]
        temp_2 <- cbind(temp[,1]-1,temp[,2:3])
        temp_3 <- cbind(temp[,1]+1,temp[,2:3])
        colnames(temp_2) <- colnames(temp)
        colnames(temp_3) <- colnames(temp)
        temp <- rbind(temp_2,temp,temp_3)
        out <- smooth.spline(temp[,1],temp[,2],control.spar=list(low=0.2,high=1))
        temp <- predict(out,(1:N)/N)$y
        sm[ii,jj,] <- (temp - min(temp)) / (max(temp) - min(temp))
        amps[ii,jj] <- max(temp) - min(temp)
    }
}


## get rid of any suspiciously large amplitudes
to_use <- rowSums(amps > 1.5) < 0.5
amps <- amps[to_use,]
sm <- sm[to_use,,]
periods <- periods[to_use]
lcs <- lcs[to_use]

amp <- median(amps)/2




######## 2. TIMES AND ERRORS FROM PANSTARRS DATA
## read data
f <- list.files("PS1_sample_LCs",full.names=TRUE)
tms <- list()
for(ii in 1:length(f)){
    temp <- read.table(f[ii],header=TRUE)
    tms[[ii]] <- data.frame(temp[,3],0,temp[,1])
    names(tms[[ii]]) <- c("time","mag","error")
}
summary(vapply(tms,nrow,0))



gammaf <- function(t,cc=0.75){
    t <- t %% 1
    return(((-2/cc)*t + 1)*(t < cc) + ((2/(1-cc))*t - (1+cc)/(1-cc))*(t >= cc))
}

## simulate magnitudes
periods <- sample(stars[,3],size=length(tms),replace=TRUE)
for(ii in 1:length(tms)){
    ## phi <- runif(1,min=0,max=2*pi)
    ## tms[[ii]][,2] <- amp*sin(((2*pi)/periods[ii]) * tms[[ii]][,1] + phi) + rnorm(nrow(tms[[ii]]),mean=0,sd=tms[[ii]][,3])
    phi <- runif(1,min=0,max=1)
    tms[[ii]][,2] <- amp*gammaf((tms[[ii]][,1]/periods[ii]) + phi) + rnorm(nrow(tms[[ii]]),mean=0,sd=tms[[ii]][,3])
    tms[[ii]][,3] <- 1
}

######## 3. WRITE OUT RESULTS
save(tms,periods,file="tms.RData")
