rm(list=ls())
#unlink("*.RData")
source('../common/funcs.R')
source('plot_options.R')
source('../common/plot_funcs.R')
source('../fit_template/fit_template.R')
load("../data/clean/sdss_rrab.RData")
library(RColorBrewer)
library(zoo)
library(ellipse)
library(xtable)

plot_foldername <- "figs"

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




## ## compute mean mag in each band / lc
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


## absolute mag period dependence
betasM <- matrix(0,nrow=length(periods),ncol=5)
for(ii in 1:length(periods)){
    betasM[ii,] <- rrmag$c0 + rrmag$p1*(log10(periods[ii]) + 0.2) + rrmag$p2*(log10(periods[ii]) + 0.2)^2
}
colnames(betasM) <- rrmag$bnd


## m <- matrix(0,nrow=length(tms),ncol=length(bands))
## for(ii in 1:length(bands)){
##     m[,ii] <- vapply(tms,function(x){mean(x[[bands[ii]]][,2])},c(0))
## }





####### BEGIN co UPDATE

## estimate ebv and mu for each lc
ests <- data.frame(ID=names(tms),mu_est=0,extcr_est=0)

## compare dust to sesar dust
m_shift <- m - betasM
rss <- rep(0,nrow(m))
for(ii in 1:nrow(m)){
    lm_fit <- lm(m_shift[ii,]~dust)
    ests[ii,2:3] <- lm_fit$coefficients
    rss[ii] <- sum(lm_fit$residuals^2)/(length(dust)-2)
}
ests[,3] <- ests[,3]*dust['r']



kpc_to_mu <- function(kpc) 5*(log10(kpc*1000)-1)
mu_to_kpc <- function(mu) (10^(mu/5 + 1)) / 1000

tab <- read.table("../data/raw/apj326724t3_mrt.txt",skip=30)
tab <- tab[,c(1,4,5)]
tab[,3] <- kpc_to_mu(tab[,3])
names(tab) <- c("ID","extcr","mu")
tab$ID <- paste0("LC_",tab$ID,".dat")

comp <- merge(tab,ests)
plot(comp$mu,comp$mu_est)
plot(comp$extcr,comp$extcr_est)


pdf("extcr_versus_mu.pdf")
y <- comp$mu_est - comp$mu
x <- comp$extcr_est-comp$extcr
par(mar=c(5,5,1,1))
plot(x,y,xlab="model estimate r Extinction - Schlegel r Extinction",ylab="model estimate mu - Sesar mu",cex.lab=1.3)
lm_fit <- lm(y ~ x)
abline(lm_fit$coeff)
abline(v=0)
abline(h=0)
points(mean(x),mean(y),pch=19,col='blue',cex=1.5)
X <- cbind(dust/dust['r'],1)
covx <- solve(t(X)%*%X)
covx <- (var(x)/covx[1,1])*covx
points(ellipse(covx,centre=c(mean(x),mean(y))),col='red',type='l',pch=2,lwd=2)
dev.off()

print(paste0("The mean error in r-band extinction is :",mean(x)))
print(paste0("The mean error in mu is :",mean(y)))




tab_out <- cbind(rrmag[,c("c0","sigc0")],0,0)
rownames(tab_out) <- rrmag$bnd
colnames(tab_out) <- c("Beta0 Original","Beta0 Uncertainty","Beta0 New","Number s.d.")

## update c0
c0_old <- rrmag$c0
cc <- mean(x) / dust['r']
rrmag$c0 <- rrmag$c0 + cc*dust


tab_out[,3] <- rrmag$c0
tab_out[,4] <- abs((tab_out[,1] - tab_out[,3]) / tab_out[,2])
tab_out <- round(tab_out,3)


tab_out <- t(tab_out)
tab_out <- tab_out[,c("u","g","r","i","z")]
out <- print(xtable(tab_out,digits=3),hline.after=0)



####### end c0 update



### REMOVE MEAN AND PHASE ALIGN LIGHTCURVES

## make lc_grid mean 0 for each band, source
for(ii in 1:Nlc){
    for(jj in 1:5){
        lc_grid[ii,,jj] <- lc_grid[ii,,jj] - mean(lc_grid[ii,,jj])
    }
}


## plot 1 function
ix_high <- 2
lc <- TMtoLC(tms[[ix_high]])
lc[,1] <- (lc[,1] %% periods[ix_high]) / periods[ix_high]
pdf(paste0(plot_foldername,"/lc_orig.pdf"),height=5,width=10)
par(mar=c(5,5,1,1))
plot(0,0,col=0,xlim=c(0,1),ylim=rev(range(lc[,3])),xaxs='i',
     xlab="Phase",ylab="Magnitude",cex.lab=1.5,cex.axis=1.5)
bands_order <- c("u","g","r","i","z")
for(ii in bands_order){
    lc_temp <- lc[lc$band==ii,]
    points(lc_temp[,1],lc_temp[,3],col=bandcol[ii],lwd=3)
    segments(lc_temp[,1],lc_temp[,3]+lc_temp[,4],lc_temp[,1],lc_temp[,3]-lc_temp[,4],col='grey')
}
legend("bottomleft",bands_order,col=bandcol[bands_order],cex=1.5,pch=1)
dev.off()
pdf(paste0(plot_foldername,"/lc_new.pdf"),height=5,width=10)
par(mar=c(5,5,1,1))
plot(0,0,col=0,xlim=range(t),ylim=rev(range(lc_grid[ix_high,,])),xaxs='i',
     xlab="Phase",ylab="Normalized Magnitude",cex.lab=1.5,cex.axis=1.5)
bands_order <- c("u","g","r","i","z")
for(ii in bands_order){
    points(t,lc_grid[ix_high,,ii],type='l',col=bandcol[ii],lwd=3)
}
legend("bottomleft",bands_order,col=bandcol[bands_order],lwd=4,cex=1.5)
dev.off()



bands_order <- c("u","g","r","i","z")
ix_high <- 2
t <- (1:100)/100
for(jj in bands_order){
    ylim <- quantile(lc_grid[,,jj],c(.999,.001))
    pdf(paste0(plot_foldername,"/unaligned_",jj,".pdf"),height=5,width=10)
    par(mar=c(5,5,1,1))
    plot(0,0,xlim=range(t),ylim=ylim,col=0,xlab="Phase",ylab=paste0(jj," Band Mag"),xaxs='i',cex.lab=1.9,cex.axis=1.9)
    for(ii in 1:Nlc){
        points(t,lc_grid[ii,,jj],type='l',col="#00000030")
    }
    points(t,lc_grid[ix_high,,jj],type='l',col='red',lwd=4)
    dev.off()
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


bands_order <- c("u","g","r","i","z")
ix_high <- 2
t <- (1:100)/100
for(jj in bands_order){
    ylim <- quantile(lc_grid[,,jj],c(.999,.001))
    pdf(paste0(plot_foldername,"/aligned_grid",jj,".pdf"),height=5,width=10)
    par(mar=c(5,5,1,1))
    plot(0,0,xlim=range(t),ylim=ylim,col=0,xlab="Phase",ylab=paste0(jj," Band Mag"),xaxs='i',cex.lab=1.9,cex.axis=1.9)
    for(ii in 1:Nlc){
        points(t,lc_grid[ii,,jj],type='l',col="#00000030")
    }
    points(t,lc_grid[ix_high,,jj],type='l',col='red',lwd=4)
    dev.off()
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



bands_order <- c("u","g","r","i","z")
ix_high <- 2
t <- (1:100)/100
for(jj in bands_order){
    ylim <- quantile(lc_grid[,,jj]/out$a,c(.999,.001))
    pdf(paste0(plot_foldername,"/aligned_amp",jj,".pdf"),height=5,width=10)
    par(mar=c(5,5,1,1))
    plot(0,0,xlim=range(t),ylim=ylim,col=0,xlab="Phase",ylab=paste0(jj," Band Mag"),xaxs='i',cex.lab=1.9,cex.axis=1.9)
    for(ii in 1:Nlc){
        points(t,lc_grid[ii,,jj]/out$a[ii],type='l',col="#00000030")
    }
    points(t,lc_grid[ix_high,,jj]/out$a[ix_high],type='l',col='red',lwd=4)
    dev.off()
}



bands_order <- c("u","g","r","i","z")
t <- (1:100)/100
for(jj in bands_order){
    ylim <- quantile(lc_grid[,,jj]/out$a,c(.999,.001))
    pdf(paste0(plot_foldername,"/aligned_amp_median",jj,".pdf"),height=5,width=10)
    par(mar=c(5,5,1,1))
    plot(0,0,xlim=range(t),ylim=ylim,col=0,xlab="Phase",ylab=paste0(jj," Band Mag"),xaxs='i',cex.lab=1.9,cex.axis=1.9)
    for(ii in 1:Nlc){
        points(t,lc_grid[ii,,jj]/out$a[ii],type='l',col="#00000030")
    }
    points(t,out$Y[,jj],type='l',col=bandcol[jj],lwd=4)
    dev.off()
}




## reorder output so templates is a B x T matrix
templates <- t(out$Y)
rownames(templates) <- dimnames(lc_grid)[[3]]
## rescale so g-band peak-to-peak amplitude is 1
gscale <- diff(range(templates["g",]))
templates <- templates / gscale
amps <- out$a*gscale
## phase shift so g-band max is phase 0
ix <- which.min(templates['g',])
if(ix != 1){
    templates <- cbind(templates[,ix:ncol(templates)],templates[,1:(ix-1)])
}

ComputeDerivative <- function(x,len,gap=2){
    n <- length(x)
    x <- c(x[(n-gap+1):n],x,x[1:gap])
    return((x[(2*gap+1):(n+2*gap)] - x[1:n]) / (2*gap*len))
}

len <- t[2] - t[1]
templatesd <- t(apply(templates,1,function(x){ComputeDerivative(x,len)}))




tem <- list(dust=dust,
            templates=templates,
            templatesd=templatesd)

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



## #### TEMPLATE MODEL refinement
## #### fit model on all light curves, compute residuals, refit templates, recompute templatesd
## #### this addresses oversmoothing at local min/max and produces better fits

## coeffs <- matrix(0,nrow=length(tms),ncol=4)
## colnames(coeffs) <- c("mu","ebv","amp","phase")
## for(ii in 1:length(tms)){
##     lc <- TMtoLC(tms[[ii]])
##     p_est <- periods[ii]
##     omega_est <- 1/p_est
##     coeffs[ii,] <- ComputeCoeffs(lc,omega_est,tem,use.errors=FALSE)
## }



## lcs <- lapply(tms,TMtoLC)
## lcs_resid <- lcs
## for(ii in 1:length(lcs)){
##     omega_est <- 1/periods[ii]
##     lcs_resid[[ii]][,1] <- (lcs[[ii]][,1]*omega_est + coeffs[ii,4]) %% 1.0
##     lcs_resid[[ii]][,3] <- lcs[[ii]][,3] - PredictTimeBand(lcs[[ii]][,1],lcs[[ii]][,2],omega_est,coeffs[ii,],tem)
## }

## lcs_resid <- do.call(rbind,lcs_resid)

## bs <- unique(lcs_resid$band)

## abs_mag_shift <- rep(0,length(bs))
## names(abs_mag_shift) <- names(bs)




## for(ii in 1:length(bs)){
##     lcs_resid_band <- lcs_resid[lcs_resid$band==bs[ii],]
##     abs_mag_shift[ii] <- mean(lcs_resid_band[,3])
##     pdf(paste0(plot_foldername,"/sdss_iter1_residuals_",bs[ii],".pdf"))
##     plot(0,0,ylim=c(.3,-.3),
##          xlab="phase",ylab="magnitude residual",
##          xlim=c(0,2),xaxs='i',col=0,main=paste0(bs[ii]," band residuals"))
##     lc1 <- lcs_resid_band
##     lc1 <- lc1[order(lc1[,1]),]
##     lc2 <- lc1
##     lc2[,1] <- lc1[,1] + 1
##     lc_temp <-rbind(lc1,lc2)
##     points(lc_temp$time,lc_temp$mag,
##            col="#00000030")

##     out <- rollmean(lc_temp[,3],201,na.pad=TRUE)
##     abline(h=0,lwd=2)
##     abline(h=abs_mag_shift[ii],col='blue',lwd=2)
##     points(lc_temp[,1],out,type='l',col='red',lwd=2)
##     dev.off()
## }

## library(splines)


## nknots <- 200
## knots <- seq(0,2,length.out=nknots+2)[2:(nknots+1)]
## coeffs.spline <- matrix(0,nrow=length(lc_temp),ncol=nknots+2)
## preds <- ns(lc_temp[,1],knots=knots,Boundary.knots=c(0,2))
## lm.fit <- lm(lc_temp[,3]~preds)
## spline_coeffs <- lm.fit$coefficients
## plot(lc_temp[,1],lc_temp[,3],xlab="phase",ylab="brightness (mag)",col="#00000010",ylim=c(.3,-.3))
## points(lc_temp[,1],matrix(c(rep(1,nrow(preds)),preds),nrow=nrow(preds))%*%spline_coeffs,
##        type='l',lwd=2,col='red')


## out <- matrix(c(rep(1,nrow(preds)),preds),nrow=nrow(preds))%*%spline_coeffs

## #### TODO: from lc1, estimate deviation from template as a function of phase
## decLocations <- tem$temp_time[1:(length(tem$temp_time)-1)]
## dec <- findInterval(lc1[,1],c(decLocations, 1))
## out <- tapply(lc1[,3],INDEX=as.factor(dec),FUN=mean)
## out[length(out)+1] <- out[1]
## plot(out) ## looks very bad

## plot(tem$templates[1,])
## a

### should we scale the template correction by a?

## returns matrix with nrow=length(p) and ncol=ncol(tem$betas)
## rows correspond to periods in p and columns are absolute
## magnitudes
tem$betas <- t(rrmag[,c("c0","p1","p2")])
colnames(tem$betas) <- rrmag$bnd
tem$abs_mag <- function(p,tem){
    X <- cbind(1,log10(p) + 0.2,(log10(p) + 0.2)^2)
    return(X%*%tem$betas)
}


## set model error initially to 0, so subsequent code runs
tem$model_error <- rep(0,length(tem$dust))
names(tem$model_error) <- names(tem$dust)


## TODO: improve this by using different level of error
## in each filter, need to trust errors more before doing this
med_res <- matrix(0,ncol=length(bands),nrow=length(tms))
for(ii in 1:nrow(med_res)){
    lc <- TMtoLC(tms[[ii]])
    coeffs <- ComputeCoeffs(lc,1/periods[ii],tem,use.errors=FALSE)
    preds <- PredictTimeBand(lc[,1],lc[,2],1/periods[ii],coeffs,tem)
    a <- tapply((preds - lc[,3])^2 - lc[,4]^2,INDEX=lc[,2],FUN=mean)
    med_res[ii,] <- a
}
model_error <- sqrt(apply(med_res,2,mean))
names(model_error) <- bands
tem$model_error <- model_error
tem$model_error[] <- mean(tem$model_error) ## makes model error same in all filters



### VISUALIZE TEMPLATES
## templates only
## colors for plotting
ylim <- c(-.7,.4)
xlim <- range(c(t,t+1))
bands_order <- c("u","g","r","i","z")
pdf(paste0(plot_foldername,"/sdss_templates.pdf"),height=8,width=15)
par(mar=c(5,5,1,1))
plot(0,0,col=0,ylim=rev(ylim),xlim=xlim,xlab="Phase",ylab=expression("Normalized Mag"~gamma),cex.lab=1.5,xaxs='i')
for(ii in bands_order){
    points(c(t,t+1),c(templates[ii,],templates[ii,]),type='l',lwd=4,
           col=bandcol[ii])
}
legend("bottomleft",bands_order,col=bandcol[bands_order],lwd=4,cex=1.5)
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





### compute residuals as a function of band
lcs <- lapply(tms,TMtoLC)
lcs_resid <- lcs
for(ii in 1:length(lcs)){
    omega_est <- 1/periods[ii]
    coeffs <- ComputeCoeffs(lcs[[ii]],omega_est,tem)
    lcs_resid[[ii]][,1] <- (lcs[[ii]][,1]*omega_est + coeffs[4]) %% 1.0
    lcs_resid[[ii]][,3] <- (lcs[[ii]][,3] - PredictTimeBand(lcs[[ii]][,1],lcs[[ii]][,2],omega_est,coeffs,tem) +
                            tem$abs_mag(1/omega_est,tem)[1,][lcs[[ii]]$band])
}

tms_resid <- lapply(lcs_resid,LCtoTM)



bands <- unique(unlist(lapply(tms_resid,names)))
betas_est <- matrix(0,nrow=3,ncol=length(bands))
rownames(betas_est) <- rownames(tem$betas)
colnames(betas_est) <- bands

pdf("period_mag_sdss.pdf",width=12,height=8)
par(mfcol=c(2,3),mar=c(5,5,3,1))
for(ii in 1:length(bands)){
    mags <- vector("numeric",length(tms_resid))
    band <- bands[ii]
    for(jj in 1:length(mags)){
        temp <- tms_resid[[jj]][[band]]
        if(is.null(temp)){
            mags[jj] <- NA
        } else {
            mags[jj] <- median(temp$mag)
        }
    }
    plot(periods,mags,main=band,cex.lab=1.3,cex.main=1.5,
         ylab="Absolute Magnitude",xlab="Period")
    ## plot fits to sdss
    vec <- tem$betas[,band]
    ti <- seq(min(periods),max(periods),length.out=100)
    points(ti,vec[1] + vec[2]*log10(ti+.2) + vec[3]*log10(ti+.2)^2,type='l',col="red",lwd=2)
    ## estimate fits to sdss, plot
    X <- cbind(1,log10(periods+.2),log10(periods+.2)^2)
    betas_est[,ii] <- lm(mags~X-1)$coefficients
    points(ti,(cbind(1,log10(ti+.2),log10(ti+.2)^2)%*%betas_est[,ii,drop=FALSE])[,1],type='l',col='black',lwd=2)
}
plot(0,0,col=0,axes=F,xlab="",ylab="")
legend("center",c("SDSS Median Mag, Dust and Distance Corrected","Given Relationships for SDSS","Fits to SDSS"),
       col=c("black","red","black"),lty=c(0,1,1),pch=c(1,-1,-1),lwd=2,cex=1.3)
dev.off()




#### CONSTRUCT OLD TEMPLATE, NO PERIOD-ABS MAG DEPENDENCE
## find betas at median period, old template style
## save template
save(tem_sdss,file="../fit_template/template_sdss.RData")



