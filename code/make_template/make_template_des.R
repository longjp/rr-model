rm(list=ls())
load("../data/clean/sdss_rrab.RData")

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


lcs <- vector("list",length(tms))
for(ii in 1:length(lcs)){
    fname <- gsub("./","",cat$filename[cat$ID==sdss_id[ii]],fixed=TRUE)
    lc <- read.table(paste0("known_rr_lcs/",fname),header=TRUE,stringsAsFactors=FALSE)
    lc <- lc[,c(1,4,2,3)]
    names(lc) <- c("time","band","mag","error")
    lcs[[ii]] <- lc
}
lcs[[1]]



## estimate coefficients using sdss data + known periods

## plot folded light curve









