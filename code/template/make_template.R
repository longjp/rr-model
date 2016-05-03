rm(list=ls())

source('func.R')

library(RColorBrewer)
library(scales)

rrlyrae <- read.table("apj326724t2_mrt.txt",skip=42)
rrlyrae <- rrlyrae[rrlyrae[,2] == "ab",]

folder <- "../rrlyrae"
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


### smooth lightcurves using supersmoother
### place on equally spaced grid
N <- 100
t <- (1:N)/N
lc_grid <- array(0,c(length(tmss),N,5),dimnames=list(NULL,NULL,bands))
for(ii in 1:length(tmss)){
    for(jj in 1:length(bands)){
        lc <- tmss[[ii]][[bands[jj]]]
        lc[,1] <- (lc[,1] %% rrlyrae[ii,3]) / rrlyrae[ii,3]
        temp <- supsmu(lc[,1],lc[,2],periodic=TRUE)
        ymean <- mean(c(temp$y[1],temp$y[length(temp$y)]))
        temp$y <- c(ymean,temp$y,ymean)
        temp$x <- c(0,temp$x,1)
        lc_grid[ii,,jj] <- approx(temp$x,temp$y,xout=t,rule=2)$y
    }
}




## compute mean mag in each band / lc
m <- apply(lc_grid,c(1,3),mean)
pdf("band_means.pdf")
pairs(m)
dev.off()

## load extinction parameters
dat <- read.table("extc.dat")
dat <- dat[dat$V1=="SDSS",c(2,3)]
dust <- dat[,2]
names(dust) <- dat[,1]
dust <- dust[order(names(dust))]


## estimate d and alpha for each lc
rs <- matrix(0,nrow=nrow(m),ncol=ncol(m))
for(ii in 1:nrow(m)){
    rs[ii,] <- lm(m[ii,]~dust)$residuals
}
## betas is the mean of the residuals (see note)
betas <- colMeans(rs)
colnames(rs) <- names(dust)
rs <- t(t(rs) - betas)

#### NOTE: CORRELATION IN RESIDUALS FOR G AND U APPEARS RELATED TO PERIOD


## construct lc with overall mean, dust, and band effects removed
for(ii in 1:length(tmss)){
    for(jj in 1:5){
        lc_grid[ii,,jj] <- lc_grid[ii,,jj] - mean(lc_grid[ii,,jj]) + rs[ii,jj]
    }
}

## determine amplitude vector by computing svd of amps

amps <- apply(lc_grid,c(1,3),function(x){mean(abs(x))})

pairs(cbind(amps,period=rrlyrae[,3]))

sv <- svd(amps)
pred <- sv$d[1]*sv$u[,1,drop=FALSE]%*%matrix(sv$v[,1],ncol=5)
##pairs(amps-pred)
amps <- abs(sv$v[,1])
names(amps) <- bands

## phase (initial) and amplitude registration
## improved phase registration using squared differences next
for(ii in 1:length(tmss)){
    temp <- lc_grid[ii,,1]
    ix <- which.min(temp)
    if(ix > 1.5){
        for(jj in 1:5){
            temp <- lc_grid[ii,,jj]
            lc_grid[ii,,jj] <- c(temp[ix:N],temp[1:(ix-1)]) / pred[ii,jj]
        }
    } else {
        for(jj in 1:5){
            temp <- lc_grid[ii,,jj]
            lc_grid[ii,,jj] <- temp / pred[ii,jj]
        }
    }
}

## phase align and compute templates for each band
del <- phase(lc_grid[,,1])
x_new <- list()
templates <- matrix(0,nrow=5,ncol=N)
rownames(templates) <- bands
for(JJ in 1:5){
    n <- nrow(lc_grid[,,JJ])
    x_new[[JJ]] <- t(vapply(1:n,function(ii){phase_shift(lc_grid[ii,,JJ],del[ii])},rep(0,N)))
    templates[JJ,] <- apply(x_new[[JJ]],2,median)
}

## visualize curves with templates
## cols <- brewer.pal(10,name="RdBu")
## decLocations <- quantile(rrlyrae[,3], probs = seq(0.1,0.9,by=0.1),type=4)
## dec <- findInterval(rrlyrae[,3],c(-Inf,decLocations, Inf))
## yedge <- 2.5
## del <- phase(lc_grid[,,1])
## for(JJ in 1:5){
##     n <- nrow(lc_grid[,,JJ])
##     dev.new()
##     par(mfcol=c(2,1))
##     ylim <- range(lc_grid[,,JJ])
##     ylim <- c(-yedge,yedge)
##     plot(0,0,ylim=ylim,xlim=c(0,1),col=0,xlab="phase",ylab="mag",xaxs="i")
##     for(ii in 1:length(tmss)){
##         points(t,lc_grid[ii,,JJ],type='l',col=cols[dec[ii]])
##     }
##     ylim <- range(x_new[[JJ]])
##     ylim <- c(-yedge,yedge)
##     plot(0,0,col=0,ylim=ylim,xlim=c(0,1),xlab="",ylab="",xaxs="i")
##     for(ii in 1:nrow(x_new[[JJ]])){
##         points(t,x_new[[JJ]][ii,],type='l',col=cols[dec[ii]])
##     }
##     points(t,templates[JJ,],lwd=3,type='l')
## }

## plot of templates
## ylim <- range(templates)
## xlim <- range(t)
## pdf("templates2.pdf",width=8,height=4)
## plot(0,0,col=0,ylim=ylim,xlim=xlim)
## for(ii in 1:5){
##     points(t,templates[ii,],type='l',lwd=2)
## }
## dev.off()


## scale templates by amps, don't save amps
for(ii in 1:length(amps)){
    templates[ii,] <- amps[ii]*templates[ii,]
}


##plot of templates
ylim <- range(templates)
xlim <- range(t)
pdf("templates.pdf",height=8,width=12)
par(mar=c(5,5,1,1))
plot(0,0,col=0,ylim=rev(ylim),xlim=xlim,xlab="Phase",ylab="Normalized Mag",cex.lab=1.5)
for(ii in 1:5){
    points(t,templates[ii,],type='l',lwd=3,col=ii,lty=ii)
}
legend("bottomleft",bands,col=1:length(bands),lty=1:length(bands),lwd=2,cex=1.5)
dev.off()

ComputeDerivative <- function(x,len,gap=2){
    n <- length(x)
    x <- c(x[(n-gap+1):n],x,x[1:gap])
    return((x[(2*gap+1):(n+2*gap)] - x[1:n]) / (2*gap*len))
}

len <- t[2] - t[1]
templatesd <- t(apply(templates,1,function(x){ComputeDerivative(x,len)}))

## ##plot of template derivatives
## dev.new()
## ylim <- range(templatesd)
## xlim <- range(t)
## plot(0,0,col=0,ylim=ylim,xlim=xlim)
## for(ii in 1:5){
##     points(t,templatesd[ii,],type='l',lwd=2)
## }

tem <- list(betas=betas,dust=dust,
            templates=templates,templatesd=templatesd)
save(tem,file="make_template.RData")
