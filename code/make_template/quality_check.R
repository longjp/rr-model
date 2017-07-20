### fits well sampled light curves to templates, calculates median residuals
### used to check that templates are working well, compare minor
### implementation differences
rm(list=ls())
library(zoo)
source('../common/funcs.R')
source('../fit_template/fit_template.R')
load("../data/clean/sdss_rrab.RData")
load("../fit_template/template_sdss.RData")


## store plots in quality_check folder
unlink("quality_check",recursive=TRUE)
dir.create("quality_check")



## find model error
med_res <- vector("numeric",length(tms))
for(ii in 1:length(med_res)){
    lc <- TMtoLC(tms[[ii]])
    coeffs <- ComputeCoeffs(lc,1/periods[ii],tem)
    preds <- PredictTimeBand(lc[,1],lc[,2],1/periods[ii],coeffs,tem)
    med_res[ii] <- median(abs(preds - lc[,3]))
}


print("total model error:")
hist(med_res)
median(med_res) ## 0.030333 seems reasonable for this number


## find model error caused by shape and fixing amplitude
med_res <- vector("numeric",length(tms))
for(ii in 1:length(med_res)){
    temp <- tms[[ii]]
    for(jj in 1:length(temp)){
        temp[[jj]][,2] <- temp[[jj]][,2] - mean(temp[[jj]][,2]) + tem$betas[jj]
    }
    lc <- TMtoLC(temp)
    coeffs <- ComputeCoeffs(lc,1/periods[ii],tem)
    preds <- PredictTimeBand(lc[,1],lc[,2],1/periods[ii],coeffs,tem)
    med_res[ii] <- median(abs(preds - lc[,3]))
}

print("model error caused by shape and fixing amplitude ratio:")
hist(med_res)
median(med_res) ## perfect mean offers only small advantage over actual model


## find model error due to shape by fitting
## model individually for each band/lc
## this removes model error caused by mean and fixing amplitude ratio
med_res <- vector("numeric",length(tms))
for(ii in 1:length(med_res)){
    lc <- TMtoLC(tms[[ii]])
    bands <- unique(lc$band)
    errors <- list()
    for(jj in 1:length(bands)){
        temp <- lc[lc$band %in% bands[jj],]
        coeffs <- ComputeCoeffs(temp,1/periods[ii],tem,use.dust=FALSE)
        preds <- PredictTimeBand(temp[,1],temp[,2],1/periods[ii],coeffs,tem)
        errors[[jj]] <- abs(preds - temp[,3])
    }
    med_res[ii] <- median(unlist(errors))
}

print("model error caused only by shape:")
hist(med_res,main="shape only")
median(med_res) 


## how much error due only to photometry?
med_res <- vector("numeric",length(tms))
for(ii in 1:length(med_res)){
    lc <- TMtoLC(tms[[ii]])
    med_res[ii] <- median(abs(rnorm(nrow(lc),mean=0,sd=lc[,4])))
}

print("photometric error:")
hist(med_res)
median(med_res)






## u band has .044 error, the worst

## are templates mean 0?
print("the mean of each template is:")
rowMeans(tem$templates)





### TODO: clean below here
######### analyze residuals as a function of phase (should be mean 0)

## compute residuals
lcs <- lapply(tms,TMtoLC)
coeffs <- matrix(0,nrow=length(tms),ncol=4)
colnames(coeffs) <- c("mu","ebv","amp","phase")
for(ii in 1:length(tms)){
    p_est <- periods[ii]
    omega_est <- 1/p_est
    coeffs[ii,] <- ComputeCoeffs(lcs[[ii]],omega_est,tem)
}


lcs_resid <- lcs
for(ii in 1:length(lcs)){
    omega_est <- 1/periods[ii]
    lcs_resid[[ii]][,1] <- (lcs[[ii]][,1]*omega_est + coeffs[ii,4]) %% 1.0
    lcs_resid[[ii]][,3] <- lcs[[ii]][,3] - PredictTimeBand(lcs[[ii]][,1],lcs[[ii]][,2],omega_est,coeffs[ii,],tem)
}



lcs_resid <- do.call(rbind,lcs_resid)



for(ii in 1:length(bands)){
    lcs_resid_band <- lcs_resid[lcs_resid$band==bands[ii],]
    pdf(paste0("quality_check/residuals_",bands[ii],".pdf"))
    plot(0,0,ylim=c(.3,-.3),
         xlab="phase",ylab="magnitude residual",
         xlim=c(0,2),xaxs='i',col=0,main=paste0(bands[ii]," band residuals"))
    lc1 <- lcs_resid_band
    lc1 <- lc1[order(lc1[,1]),]
    lc2 <- lc1
    lc2[,1] <- lc1[,1] + 1
    lc_temp <-rbind(lc1,lc2)
    points(lc_temp$time,lc_temp$mag,
           col="#00000030")

    out <- rollmedian(lc_temp[,3],51,na.pad=TRUE)
    abline(h=0,lwd=2)
    points(lc_temp[,1],out,type='l',col='red',lwd=2)
    dev.off()
}

## segments(lc_temp$time,
    ##          lc_temp$mag+lc_temp$error,
    ##          lc_temp$time,
    ##          lc_temp$mag-lc_temp$error)

