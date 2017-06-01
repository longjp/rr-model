rm(list=ls())
unlink("*.RData")
unlink("*.pdf")
source('../common/funcs.R')
load("../data/clean/sdss_rrab.RData")
library(RColorBrewer)

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

### smooth lightcurves using supersmoother
### place on equally spaced grid
N <- 100
t <- (1:N)/N
lc_grid <- array(0,c(Nlc,N,5),dimnames=list(NULL,NULL,bands))
for(ii in 1:Nlc){
    for(jj in 1:length(bands)){
        lc <- tms[[ii]][[bands[jj]]]
        lc[,1] <- (lc[,1] %% periods[ii]) / periods[ii]
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

## load and compute beta (M_x) parameters using median of periods
rrmag <- read.table("rrmag.dat",header=TRUE)
rrmag <- rrmag[rrmag$Sys=="SDSS",]
rrmag <- rrmag[order(rrmag$bnd),]
lpmed <- log10(median(periods))
betas <- rrmag$c0 + rrmag$p1*(lpmed + 0.2) + rrmag$p2*(lpmed + 0.2)^2
names(betas) <- rrmag$bnd

## estimate ebv and mu for each lc
rs <- matrix(0,nrow=nrow(m),ncol=ncol(m))
m_shift <- t(t(m) - betas)
for(ii in 1:nrow(m)){
    rs[ii,] <- lm(m_shift[ii,]~dust)$residuals
}

#### NOTE: CORRELATION IN RESIDUALS FOR G AND U APPEARS RELATED TO PERIOD

## construct lc with overall mean, dust, and band effects removed
for(ii in 1:Nlc){
    for(jj in 1:5){
        lc_grid[ii,,jj] <- lc_grid[ii,,jj] - mean(lc_grid[ii,,jj]) ## + rs[ii,jj]
    }
}

## determine amplitude vector by computing svd of amps
amps <- apply(lc_grid,c(1,3),function(x){mean(abs(x))})

##pairs(cbind(amps,periods))
sv <- svd(amps)
pred <- sv$d[1]*sv$u[,1,drop=FALSE]%*%matrix(sv$v[,1],ncol=5)
##pairs(amps-pred)
amps <- abs(sv$v[,1])
names(amps) <- bands


### first phase align light curves
## crude phase aligns, make minimums = phase 0
for(ii in 1:Nlc){
    temp <- lc_grid[ii,,1]
    ix <- which.min(temp)
    if(ix > 1.5){
        for(jj in 1:5){
            temp <- lc_grid[ii,,jj]
            lc_grid[ii,,jj] <- c(temp[ix:N],temp[1:(ix-1)])
        }
    } 
}
## sophisticated phase align, based on FDA book, see phase, phase_shift functions
del <- phase(lc_grid[,,1])
x_new <- list()
for(JJ in 1:5){
    n <- nrow(lc_grid[,,JJ])
    x_new[[JJ]] <- t(vapply(1:n,function(ii){phase_shift(lc_grid[ii,,JJ],del[ii])},rep(0,N)))
}


X <- lc_grid
for(ii in 1:dim(lc_grid)[1]){
    for(jj in 1:dim(lc_grid)[3]){
        X[ii,,jj] <- phase_shift(lc_grid[ii,,jj],del[ii])
    }
}


## TODO: error check this
SolveAGamma <- function(X,a=NULL,N=1000){
    n <- dim(X)[1]
    ## initialize a as normalized unit vector
    if(is.null(a)){
        a <- rep(1,n)
        a <- a / sqrt(sum(a*a))
    }
    ## update a and Y iteratively, N times
    for(jj in 1:N){
        Y <- vapply(1:n,function(ix){a[ix]*X[ix,,]},matrix(0,nrow=dim(X)[2],ncol=dim(X)[3]))
        Y <- apply(Y,1:2,sum)
        a <- vapply(1:n,function(ix){sum(Y*X[ix,,])},c(0))
        a <- a / sqrt(sum(a*a))
    }
    return(list(a=a,Y=Y))
}
out <- SolveAGamma(X)

Y <- out$Y


##### LOOKS BEAUTIFUL!
ylim <- range(Y)
plot(0,0,xlim=range(t),ylim=ylim,col=0)
for(ii in 1:dim(Y)[2]){
    points(t,Y[,ii],type='l')
}


## TODO: replace next two steps with optimization algorithm

## amplitude align, based on svd
for(jj in 1:5){
    for(ii in 1:nrow(pred)){
        x_new[[jj]][ii,] <- x_new[[jj]][ii,] / pred[ii,jj]
    }
}
## make templates
templates <- matrix(0,nrow=5,ncol=N)
rownames(templates) <- bands
for(JJ in 1:5){
    templates[JJ,] <- apply(x_new[[JJ]],2,median)
}




## visualize curves with templates
cols <- brewer.pal(10,name="RdBu")
decLocations <- quantile(periods, probs = seq(0.1,0.9,by=0.1),type=4)
dec <- findInterval(periods,c(-Inf,decLocations, Inf))
yedge <- 2.5
del <- phase(lc_grid[,,1])
for(JJ in 1:5){
    n <- nrow(lc_grid[,,JJ])
    dev.new()
    par(mfcol=c(2,1))
    ylim <- range(lc_grid[,,JJ])
    ylim <- c(-yedge,yedge)
    plot(0,0,ylim=ylim,xlim=c(0,1),col=0,xlab="phase",ylab="mag",xaxs="i")
    for(ii in 1:Nlc){
        points(t,lc_grid[ii,,JJ],type='l',col=cols[dec[ii]])
    }
    ylim <- range(x_new[[JJ]])
    ylim <- c(-yedge,yedge)
    plot(0,0,col=0,ylim=ylim,xlim=c(0,1),xlab="",ylab="",xaxs="i")
    for(ii in 1:nrow(x_new[[JJ]])){
        points(t,x_new[[JJ]][ii,],type='l',col=cols[dec[ii]])
    }
    points(t,templates[JJ,],lwd=3,type='l')
}

##plot of templates
ylim <- range(templates)
xlim <- range(t)
pdf("templates2.pdf",width=8,height=4)
plot(0,0,col=0,ylim=ylim,xlim=xlim)
for(ii in 1:5){
    points(t,templates[ii,],type='l',lwd=2)
}
dev.off()


## scale templates by amps, don't save amps
for(ii in 1:length(amps)){
    templates[ii,] <- amps[ii]*templates[ii,]
}

## make g band peak-to-peak amplitude 1
gscale <- diff(range(templates["g",]))
for(ii in 1:length(amps)){
    templates[ii,] <- templates[ii,] / gscale
}



##plot of templates
ylim <- range(templates)
xlim <- range(t)
pdf("templates.pdf",height=8,width=12)
par(mar=c(5,5,1,1))
plot(0,0,col=0,ylim=rev(ylim),xlim=xlim,xlab="Phase",ylab=expression("Normalized Mag"~gamma),cex.lab=1.5,xaxs='i')
for(ii in 1:5){
    points(t,templates[ii,],type='l',lwd=4,col=ii,lty=ii)
}
legend("bottomleft",bands,col=1:length(bands),lty=1:length(bands),lwd=4,cex=1.5)
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

## functions for interpolating templates
temp_time <- seq(0,1,length.out=ncol(tem$templates))
tem$temp_time <- temp_time
tem$template_funcs <- list()
for(jj in 1:nrow(tem$templates)){
    tem$template_funcs[[jj]] <- approxfun(temp_time,tem$templates[jj,])
}
names(tem$template_funcs) <- bands
      
tem$templated_funcs <- list()
for(jj in 1:nrow(tem$templatesd)){
    tem$templated_funcs[[jj]] <- approxfun(temp_time,tem$templatesd[jj,])
}
names(tem$templated_funcs) <- bands

save(tem,file="../fit_template/template.RData")
