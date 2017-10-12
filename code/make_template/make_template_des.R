## modifies sdss template to create des template
rm(list=ls())
load("../data/clean/des.RData")
load("../fit_template/template_sdss.RData")
source("../fit_template/fit_template.R")
source("../common/funcs.R")
source("../common/plot_funcs.R")
source("funcs.R")
library(zoo)
library(L1pack)

## for storing plots
plot_foldername <- "figs"


lcs_sdss <- lapply(tms_sdss,TMtoLC)
lcs_des <- lapply(tms_des,TMtoLC)

## estimate coefficients using sdss data + known periods
coeffs <- matrix(0,nrow=length(lcs_sdss),ncol=4)
colnames(coeffs) <- c("mu","ebv","amp","phase")
for(ii in 1:length(tms_sdss)){
    lc <- lcs_sdss[[ii]]
    p_est <- periods[ii]
    omega_est <- 1/p_est
    coeffs[ii,] <- ComputeCoeffs(lc,omega_est,tem)
}


##### create DES template
#####
### add Y band template that is identical to z
tem <- AddBand(tem,"Y","z")


## plot folded light curve for des, sdss
ii <- 6
par(mfcol=c(2,1))
plotLC(lcs_sdss[[ii]],periods[ii],coeffs[ii,],tem)
plotLC(lcs_des[[ii]],periods[ii],coeffs[ii,],tem)




## set dust to DES values
extc <- read.table("extc.dat",stringsAsFactors=FALSE)
extc <- extc[extc[,1]=="DES",c(2,3)]
tem$dust[extc$V2] <- extc$V3






###
### DETERMINE PERIOD-ABS MAG DEPENDENCE FOR DES BANDS
###

### compute residuals as a function of band
lcs_resid <- lcs_des
for(ii in 1:length(lcs_des)){
    omega_est <- 1/periods[ii]
    lcs_resid[[ii]][,1] <- (lcs_des[[ii]][,1]*omega_est + coeffs[ii,4]) %% 1.0
    lcs_resid[[ii]][,3] <- (lcs_des[[ii]][,3] - PredictTimeBand(lcs_des[[ii]][,1],lcs_des[[ii]][,2],omega_est,coeffs[ii,],tem) +
                            tem$abs_mag(1/omega_est,tem)[1,][lcs_des[[ii]]$band])
}

tms_resid <- lapply(lcs_resid,LCtoTM)

### plot residuals as a function of band and estimate function
bands <- unique(unlist(lapply(tms_resid,names)))
betas_des <- matrix(0,nrow=3,ncol=length(bands))
rownames(betas_des) <- rownames(tem$betas)
colnames(betas_des) <- bands

pdf("period_mag_des.pdf",width=12,height=8)
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
    plot(periods,mags,main=band,ylim=c(.25,.9),cex.lab=1.3,cex.main=1.5,
         ylab="Absolute Magnitude",xlab="Period")
    ## plot fits to sdss
    vec <- tem$betas[,band]
    ti <- seq(min(periods),max(periods),length.out=100)
    points(ti,vec[1] + vec[2]*(log10(ti)+.2) + vec[3]*(log10(ti)+.2)^2,type='l',col="red",lwd=2)
    ## estimate fits to des, plot
    X <- cbind(1,log10(periods)+.2,(log10(periods)+.2)^2)
    betas_des[,ii] <- lm(mags~X-1)$coefficients
    points(ti,(cbind(1,log10(ti)+.2,(log10(ti)+.2)^2)%*%betas_des[,ii,drop=FALSE])[,1],type='l',col='black',lwd=2)
}
plot(0,0,col=0,axes=F,xlab="",ylab="")
legend("center",c("DES Median Mag, Dust and Distance Corrected","Given Relationships for SDSS","Fits to DES"),
       col=c("black","red","black"),lty=c(0,1,1),pch=c(1,-1,-1),lwd=2,cex=1.3)
dev.off()





## store new DES values for betas
tem$betas <- betas_des



## lcs_resid <- do.call(rbind,lcs_resid)

## bs <- unique(lcs_resid$band)

## abs_mag_shift <- rep(0,length(bs))
## names(abs_mag_shift) <- names(bs)

## for(ii in 1:length(bs)){
##     lcs_resid_band <- lcs_resid[lcs_resid$band==bs[ii],]
##     abs_mag_shift[ii] <- mean(lcs_resid_band[,3])
##     pdf(paste0(plot_foldername,"/residuals_",bs[ii],"_old.pdf"))
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

##     out <- rollmedian(lc_temp[,3],51,na.pad=TRUE)
##     abline(h=0,lwd=2)
##     abline(h=abs_mag_shift[ii],col='blue',lwd=2)
##     points(lc_temp[,1],out,type='l',col='red',lwd=2)
##     dev.off()
## }



## ### STEP 2: shift absolute magnitudes to make mean 0
## ## shift absolute magnitudes to mean 0
## tem$betas[bs] <- tem$betas[bs] + abs_mag_shift



## ### STEP 3: rerun step 1, residuals should now have mean 0
## lcs_resid <- lcs_des
## for(ii in 1:length(lcs_des)){
##     omega_est <- 1/periods[ii]
##     lcs_resid[[ii]][,1] <- (lcs_des[[ii]][,1]*omega_est + coeffs[ii,4]) %% 1.0
##     lcs_resid[[ii]][,3] <- lcs_des[[ii]][,3] - PredictTimeBand(lcs_des[[ii]][,1],lcs_des[[ii]][,2],omega_est,coeffs[ii,],tem)
## }

## lcs_resid <- do.call(rbind,lcs_resid)

## bs <- unique(lcs_resid$band)

## abs_mag_shift <- rep(0,length(bs))
## names(abs_mag_shift) <- names(bs)

## for(ii in 1:length(bs)){
##     lcs_resid_band <- lcs_resid[lcs_resid$band==bs[ii],]
##     abs_mag_shift[ii] <- mean(lcs_resid_band[,3])
##     pdf(paste0(plot_foldername,"/residuals_",bs[ii],"_new.pdf"))
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

##     out <- rollmedian(lc_temp[,3],51,na.pad=TRUE)
##     abline(h=0,lwd=2)
##     abline(h=abs_mag_shift[ii],col='blue',lwd=2)
##     points(lc_temp[,1],out,type='l',col='red',lwd=2)
##     dev.off()
## }





####
#### SAVE OUTPUT
####

## get rid of u band
tem <- RemoveBand(tem,"u")

## save template
save(tem,file="../fit_template/template_des.RData")





