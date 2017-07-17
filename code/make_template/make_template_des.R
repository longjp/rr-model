rm(list=ls())
load("../data/clean/sdss_rrab.RData")
load("../fit_template/template.RData")
source("../fit_template/fit_template.R")
source("../common/funcs.R")

##### make plot with sloan data in grey or very light, des darker
##### print out every one of these light curves, sanity check
##### estimate Y offset and shape, compare to z

## makes nice plot
plotLC <- function(lc,p_est,tem,coeffs=NULL,add=FALSE,pch=TRUE){
    colpch <- 1:length(tem$betas)
    if(!pch){
        colpch <- rep(length(tem$betas)+1,length(tem$betas))
    }
    names(colpch) <- names(tem$betas)
    lc1 <- lc
    lc1[,1] <- (lc$time %% p_est)/p_est
    lc2 <- lc1
    lc2[,1] <- lc1[,1] + 1
    lc_temp <-rbind(lc1,lc2)
    if(add==FALSE){
        plot(0,0,ylim=rev(range(lc_temp$mag)),
             xlab="time",ylab="magnitude",
             xlim=c(0,2),xaxs='i',col=0)
    }
    points(lc_temp$time,lc_temp$mag,
         col=colpch[lc_temp$band],pch=colpch[lc_temp$band])
    segments(lc_temp$time,
             lc_temp$mag+lc_temp$error,
             lc_temp$time,
             lc_temp$mag-lc_temp$error)
    if(!is.null(coeffs)){
        ti <- (1:100)/100
        ti <- c(ti,ti+1)
        m <- PredictAllBand(ti,1,coeffs,tem)
        for(ii in 1:length(tem$betas)){
            points(ti,m[,ii],type='l',col=colpch[names(tem$betas)[ii]])
        }
    }
}


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


## read in katelyn catalog data and get sdss ids in nice form
cat <- read.table("known_rr_lcs/rrab_only.tab",colClasses=c("numeric","character","character","character","numeric"),
                  header=TRUE)
cat$ID <- gsub(".0","",cat$ID,fixed=TRUE)
sdss_id <- gsub(".dat","",names(tms))
sdss_id <- gsub("LC_","",sdss_id)

## get l.c.s which have cross matched with des
to_use <- sdss_id %in% cat$ID
Nlc <- sum(to_use)
sdss_id <- sdss_id[to_use]
periods <- periods[to_use]
tms <- tms[to_use]
lcs_sdss <- lapply(tms,TMtoLC)


lcs_des <- vector("list",length(tms))
for(ii in 1:length(lcs_des)){
    fname <- gsub("./","",cat$filename[cat$ID==sdss_id[ii]],fixed=TRUE)
    lc <- read.table(paste0("known_rr_lcs/",fname),header=TRUE,stringsAsFactors=FALSE)
    lc <- lc[,c(1,4,2,3)]
    names(lc) <- c("time","band","mag","error")
    lcs_des[[ii]] <- lc
}




## estimate coefficients using sdss data + known periods
ii <- 1
lc_sdss <- lcs_sdss[[ii]]
lc_des <- lcs_des[[ii]]
p_est <- periods[ii]
omega_est <- 1/p_est
coeffs <- ComputeCoeffs(lc_sdss,omega_est,tem)
names(coeffs) <- c("mu","ebv","amp","phase")

## plot folded light curve
plotLC(lc_sdss,p_est,tem,coeffs)

plotLC(lc_des,p_est,tem,add=TRUE,pch=FALSE)



## compute median dust
coeffs <- matrix(0,nrow=length(tms),ncol=4)
colnames(coeffs) <- c("mu","ebv","amp","phase")
for(ii in 1:length(tms)){
    lc <- TMtoLC(tms[[ii]])
    p_est <- periods[ii]
    omega_est <- 1/p_est
    coeffs[ii,] <- ComputeCoeffs(lc,omega_est,tem)
}

hist(coeffs[,2])
summary(coeffs[,2])


## set dust to new values
extc <- read.table("extc.dat",stringsAsFactors=FALSE)
extc <- extc[extc[,1]=="DES",c(2,3)]
tem$dust[extc$V2] <- extc$V3



### TODO: at minimum DES templates should use DES dust from extc
###       changes to shape of templates may make sense as well



### how much of a problem is aliasing
times <- vector("list",length(lcs_des))
for(ii in 1:length(times)){
    temp <- sort(lcs_des[[ii]][,1])
    n <- length(temp)
    times[[ii]] <- temp[2:n] - temp[1:(n-1)]
}
times <- unlist(times)
hist(times %% 1,main="des",breaks=seq(0,1,.025))


times <- vector("list",length(lcs_des))
for(ii in 1:length(times)){
    temp <- sort(lcs_sdss[[ii]][,1])
    n <- length(temp)
    times[[ii]] <- temp[2:n] - temp[1:(n-1)]
}
times <- unlist(times)
dev.new()
hist(times %% 1,main="sdss",breaks=seq(0,1,.025))
### conclusion: aliasing will be problem with des, but not as much as with sdss
