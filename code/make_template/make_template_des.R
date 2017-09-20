## modifies sdss template to create des template
rm(list=ls())
load("../data/clean/des.RData")
load("../fit_template/template_sdss.RData")
source("../fit_template/fit_template.R")
source("../common/funcs.R")
source("../common/plot_funcs.R")
source("funcs.R")
library(zoo)

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


## plot folded light curve for des, sdss
ii <- 1
par(mfcol=c(2,1))
plotLC(lcs_sdss[[ii]],periods[ii],coeffs[ii,],tem)
plotLC(lcs_des[[ii]],periods[ii],coeffs[ii,],tem)



##### create DES template
#####
### add Y band template that is identical to z
tem <- AddBand(tem,"Y","z")

tem$betas

## set dust to DES values
extc <- read.table("extc.dat",stringsAsFactors=FALSE)
extc <- extc[extc[,1]=="DES",c(2,3)]
tem$dust[extc$V2] <- extc$V3






###
### MODEL RESIDUAL ANALYSIS
###

### TODO: reestimate c0, p1, and p2 terms for ALL bands? since these coefficients are for SDSS


### STEP 1: compute residuals as a function of phase by band
lcs_resid <- lcs_des
for(ii in 1:length(lcs_des)){
    omega_est <- 1/periods[ii]
    lcs_resid[[ii]][,1] <- (lcs_des[[ii]][,1]*omega_est + coeffs[ii,4]) %% 1.0
    lcs_resid[[ii]][,3] <- lcs_des[[ii]][,3] - PredictTimeBand(lcs_des[[ii]][,1],lcs_des[[ii]][,2],omega_est,coeffs[ii,],tem)
}

lcs_resid <- do.call(rbind,lcs_resid)

bs <- unique(lcs_resid$band)

abs_mag_shift <- rep(0,length(bs))
names(abs_mag_shift) <- names(bs)

for(ii in 1:length(bs)){
    lcs_resid_band <- lcs_resid[lcs_resid$band==bs[ii],]
    abs_mag_shift[ii] <- mean(lcs_resid_band[,3])
    pdf(paste0(plot_foldername,"/residuals_",bs[ii],"_old.pdf"))
    plot(0,0,ylim=c(.3,-.3),
         xlab="phase",ylab="magnitude residual",
         xlim=c(0,2),xaxs='i',col=0,main=paste0(bs[ii]," band residuals"))
    lc1 <- lcs_resid_band
    lc1 <- lc1[order(lc1[,1]),]
    lc2 <- lc1
    lc2[,1] <- lc1[,1] + 1
    lc_temp <-rbind(lc1,lc2)
    points(lc_temp$time,lc_temp$mag,
           col="#00000030")

    out <- rollmedian(lc_temp[,3],51,na.pad=TRUE)
    abline(h=0,lwd=2)
    abline(h=abs_mag_shift[ii],col='blue',lwd=2)
    points(lc_temp[,1],out,type='l',col='red',lwd=2)
    dev.off()
}



### STEP 2: shift absolute magnitudes to make mean 0
## shift absolute magnitudes to mean 0
tem$betas[bs] <- tem$betas[bs] + abs_mag_shift



### STEP 3: rerun step 1, residuals should now have mean 0
lcs_resid <- lcs_des
for(ii in 1:length(lcs_des)){
    omega_est <- 1/periods[ii]
    lcs_resid[[ii]][,1] <- (lcs_des[[ii]][,1]*omega_est + coeffs[ii,4]) %% 1.0
    lcs_resid[[ii]][,3] <- lcs_des[[ii]][,3] - PredictTimeBand(lcs_des[[ii]][,1],lcs_des[[ii]][,2],omega_est,coeffs[ii,],tem)
}

lcs_resid <- do.call(rbind,lcs_resid)

bs <- unique(lcs_resid$band)

abs_mag_shift <- rep(0,length(bs))
names(abs_mag_shift) <- names(bs)

for(ii in 1:length(bs)){
    lcs_resid_band <- lcs_resid[lcs_resid$band==bs[ii],]
    abs_mag_shift[ii] <- mean(lcs_resid_band[,3])
    pdf(paste0(plot_foldername,"/residuals_",bs[ii],"_new.pdf"))
    plot(0,0,ylim=c(.3,-.3),
         xlab="phase",ylab="magnitude residual",
         xlim=c(0,2),xaxs='i',col=0,main=paste0(bs[ii]," band residuals"))
    lc1 <- lcs_resid_band
    lc1 <- lc1[order(lc1[,1]),]
    lc2 <- lc1
    lc2[,1] <- lc1[,1] + 1
    lc_temp <-rbind(lc1,lc2)
    points(lc_temp$time,lc_temp$mag,
           col="#00000030")

    out <- rollmedian(lc_temp[,3],51,na.pad=TRUE)
    abline(h=0,lwd=2)
    abline(h=abs_mag_shift[ii],col='blue',lwd=2)
    points(lc_temp[,1],out,type='l',col='red',lwd=2)
    dev.off()
}





####
#### SAVE OUTPUT
####

## get rid of u band
tem <- RemoveBand(tem,"u")

## save template
save(tem,file="../fit_template/template_des.RData")





