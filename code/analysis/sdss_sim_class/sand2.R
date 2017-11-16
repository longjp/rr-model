### created for cdi grant proposal
### exploring shape of rss functions for sine model
rm(list=ls())
set.seed(1234)

## load necessary libraries and data source
library('multiband')
load("../../fit_template/template_sdss.RData")
source("../../fit_template/fit_template.R")
##source("../../common/funcs.R")
##source("../funcs.R")
load("../../data/clean/sdss_sim_class.RData")

## load some utlity functions
GetFreqs <-function(period_min,period_max,freq_del = 0.1/4000){
    freq_max <- 1/period_min
    freq_min <- 1/period_max
    return(seq(freq_min, freq_max, freq_del))
}

TMtoLC <- function(tm){
    nb <- vapply(tm,nrow,c(0))
    bs <- rep.int(names(tm),nb)
    lc <- do.call(rbind,tm)
    lc <- data.frame(time=lc[,1],band=bs,mag=lc[,2],error=lc[,3])
    lc[,2] <- as.character(lc[,2])
    return(lc)
}
## returns dust corrected lightcurve
##    arguments:
##          tm : lightcurve
##         ebv : dust for lightcurve tm
##         tem : templates, contains extinction law in tem$dust
DustCorrect <- function(tm,ebv,tem){
    if(mean(names(tm) %in% names(tem$dust))!=1){
        print("lightcurve has bands: ",names(tm))
        print("template only has extinction law for:",names(tem$dust))
        stop()
    }
    for(jj in 1:length(tm)){
        tm[[jj]][,2] <- tm[[jj]][,2] - tem$dust[names(tm)[jj]]*ebv
    }
    return(tm)
}


## parameters for simulation
NN <- 10 ## number of newton steps 
omegas <- GetFreqs(0.4,0.95) ## frequency grid

## create dust corrected tms
ebv <- extcr / tem$dust['r']
tmsc <- mapply(DustCorrect,tms,ebv,MoreArgs=list(tem=tem),SIMPLIFY=FALSE)
tmsc_FULL <- mapply(DustCorrect,tms_FULL,ebv,MoreArgs=list(tem=tem),SIMPLIFY=FALSE)

ii <- 1
p_grid <- rev(1/omegas)

col <- "#00000030"
wid <- 600
hei <- 350

## get light curve
tm <- tms_FULL[[ii]]
ix <- sample(nrow(tm[[1]]),20)
tm <- list(tm[[1]][ix,])
names(tm) <- "g"
lc <- TMtoLC(tm)
n <- nrow(lc)

#### SINE MODEL
## full
rss <- pgls(tm,periods=p_grid,BCD_flag=FALSE)$rss_ls
png("rss_sin_full.png",width=wid,height=hei)
par(mar=c(5,5,1,1))
plot(p_grid,rss/n,xlab="Period",ylab="RSS/n",cex.lab=1.3,cex.axis=1.3,xaxs='i',col=col)
abline(v=periods[ii],col='red',lwd=2)
abline(v=p_grid[which.min(rss)],col='orange',lwd=2)
dev.off()

#### RRL MODEL
rss <- FitTemplate(lc,omegas,tem,NN=NN,use.errors=TRUE,use.dust=FALSE)
png("rss_rr_model.png",width=wid,height=hei)
par(mar=c(5,5,1,1))
plot(p_grid,rev(rss)/n,xlab="Period",ylab="RSS/n",cex.lab=1.3,cex.axis=1.3,xaxs='i',col=col)
abline(v=periods[ii],col='red',lwd=2)
abline(v=p_grid[which.min(rev(rss))],col='orange',lwd=2)
dev.off()
