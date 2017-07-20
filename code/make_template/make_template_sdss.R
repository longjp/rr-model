rm(list=ls())
unlink("*.RData")
source('../common/funcs.R')
source('../fit_template/fit_template.R')
load("../data/clean/sdss_rrab.RData")
library(RColorBrewer)

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

## interpolate light curves on a grid
N <- 100
t <- (1:N)/N
lc_grid <- array(0,c(Nlc,N,5),dimnames=list(NULL,NULL,bands))
for(ii in 1:Nlc){
    for(jj in 1:length(bands)){
        lc <- tms[[ii]][[bands[jj]]]
        lc[,1] <- (lc[,1] %% periods[ii]) / periods[ii]
        ords <- order(lc[,1])
        x <- lc[ords,1]
        y <- lc[ords,2]
        ymean <- mean(c(y[1],y[length(y)]))
        y <- c(ymean,y,ymean)
        x <- c(0,x,1)
        y_approx <- approx(x,y,xout=t,rule=2)$y
        lc_grid[ii,,jj] <- y_approx
    }
}



## compute mean mag in each band / lc
m <- apply(lc_grid,c(1,3),mean)
pdf("figs/band_means.pdf")
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

## make lc_grid mean 0 for each band, source
for(ii in 1:Nlc){
    for(jj in 1:5){
        lc_grid[ii,,jj] <- lc_grid[ii,,jj] - mean(lc_grid[ii,,jj])
    }
}

### PHASE ALIGN
## method 1: crude phase aligns, make minimums = phase 0
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
## method 2: based on FDA book, see phase, phase_shift functions
del <- phase(lc_grid[,,1],niter=1000)
for(ii in 1:dim(lc_grid)[1]){
    for(jj in 1:dim(lc_grid)[3]){
        lc_grid[ii,,jj] <- phase_shift(lc_grid[ii,,jj],del[ii])
    }
}
## method 3: grid search to optimize phase (1-step)
del <- PhaseGridAll(lc_grid[,,1])
for(ii in 1:dim(lc_grid)[1]){
    for(jj in 1:dim(lc_grid)[3]){
        lc_grid[ii,,jj] <- phase_shift(lc_grid[ii,,jj],del[ii])
    }
}



####### phase aligned light curves look good
## JJ <- 0

## JJ <- JJ + 1
## ylim <- range(lc_grid[,,JJ])
## cols <- brewer.pal(10,name="RdBu")
## decLocations <- quantile(periods, probs = seq(0.1,0.9,by=0.1),type=4)
## dec <- findInterval(periods,c(-Inf,decLocations, Inf))
## plot(0,0,ylim=ylim,xlim=c(0,1),col=0,xlab="phase",ylab="mag",xaxs="i")
## for(ii in 1:Nlc){
##     points(t,lc_grid[ii,,JJ],type='l',col=cols[dec[ii]])
## }
## abline(h=0,col='red',lwd=3)

## determine amplitude and templates for sources
out <- SolveAGamma(lc_grid)

## reorder output so templates is a B x T matrix
templates <- t(out$Y)
rownames(templates) <- dimnames(lc_grid)[[3]]
## rescale so g-band peak-to-peak amplitude is 1
gscale <- diff(range(templates["g",]))
templates <- templates / gscale
amps <- out$a*gscale



### VISUALIZE TEMPLATES
## templates only
ylim <- range(templates)
xlim <- range(t)
pdf("figs/sdss_templates.pdf",height=8,width=12)
par(mar=c(5,5,1,1))
plot(0,0,col=0,ylim=rev(ylim),xlim=xlim,xlab="Phase",ylab=expression("Normalized Mag"~gamma),cex.lab=1.5,xaxs='i')
for(ii in 1:5){
    points(t,templates[ii,],type='l',lwd=4,col=ii,lty=ii)
}
legend("bottomleft",bands,col=1:length(bands),lty=1:length(bands),lwd=4,cex=1.5)
dev.off()

## templates with light curves
## cols <- brewer.pal(10,name="RdBu")
## decLocations <- quantile(periods, probs = seq(0.1,0.9,by=0.1),type=4)
## dec <- findInterval(periods,c(-Inf,decLocations, Inf))
## for(JJ in 1:5){
##     ylim <- range(templates[JJ,])
##     pdf(paste0("figs/sdss_template_",JJ,".pdf"),height=8,width=12)
##     plot(0,0,ylim=ylim,xlim=c(0,1),col=0,xlab="phase",ylab="mag",xaxs="i")
##     for(ii in 1:Nlc){
##         points(t,lc_grid[ii,,JJ]/amps[ii],type='l',col=cols[dec[ii]])
##     }
##     points(t,templates[JJ,],lwd=3,type='l')
##     dev.off()
## }

ComputeDerivative <- function(x,len,gap=2){
    n <- length(x)
    x <- c(x[(n-gap+1):n],x,x[1:gap])
    return((x[(2*gap+1):(n+2*gap)] - x[1:n]) / (2*gap*len))
}

len <- t[2] - t[1]
templatesd <- t(apply(templates,1,function(x){ComputeDerivative(x,len)}))

## VISUALIZE TEMPLATE DERIVATIVES
ylim <- range(templatesd)
xlim <- range(t)
pdf("figs/sdss_templatesd.pdf",height=8,width=12)
par(mar=c(5,5,1,1))
plot(0,0,col=0,ylim=rev(ylim),xlim=xlim,xlab="Phase",ylab=expression("Derivative Normalized Mag"~gamma),cex.lab=1.5,xaxs='i')
for(ii in 1:5){
    points(t,templatesd[ii,],type='l',lwd=4,col=ii,lty=ii)
}
legend("bottomleft",bands,col=1:length(bands),lty=1:length(bands),lwd=4,cex=1.5)
dev.off()





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


## set model error initially to 0, so subsequent code runs
tem$model_error <- rep(0,length(tem$betas))
names(tem$model_error) <- names(tem$betas)


##
## estimate model error by band
##
med_res <- matrix(0,ncol=length(bands),nrow=length(tms))
for(ii in 1:nrow(med_res)){
    lc <- TMtoLC(tms[[ii]])
    coeffs <- ComputeCoeffs(lc,1/periods[ii],tem,use.errors=FALSE)
    preds <- PredictTimeBand(lc[,1],lc[,2],1/periods[ii],coeffs,tem)
    a <- tapply((preds - lc[,3])^2 - lc[,4]^2,INDEX=lc[,2],FUN=median)
    med_res[ii,] <- a
}
model_error <- sqrt(apply(med_res,2,median))
names(model_error) <- bands
tem$model_error <- model_error




## save template
save(tem,file="../fit_template/template_sdss.RData")



