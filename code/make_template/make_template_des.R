rm(list=ls())
load("../data/clean/sdss_rrab.RData")
load("../fit_template/template.RData")
source("../fit_template/fit_template.R")
source("../common/funcs.R")
source("../common/plot_funcs.R")
library(zoo)

## for storing plots
plot_foldername <- "quality_check_des"
unlink(plot_foldername,recursive=TRUE)
dir.create(plot_foldername)

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

## read in des lcs
lcs_des <- vector("list",length(tms))
for(ii in 1:length(lcs_des)){
    fname <- gsub("./","",cat$filename[cat$ID==sdss_id[ii]],fixed=TRUE)
    lc <- read.table(paste0("known_rr_lcs/",fname),header=TRUE,stringsAsFactors=FALSE)
    lc <- lc[,c(1,4,2,3)]
    names(lc) <- c("time","band","mag","error")
    lcs_des[[ii]] <- lc
}

## estimate coefficients using sdss data + known periods
coeffs <- matrix(0,nrow=length(tms),ncol=4)
colnames(coeffs) <- c("mu","ebv","amp","phase")
for(ii in 1:length(tms)){
    lc <- lcs_sdss[[ii]]
    p_est <- periods[ii]
    omega_est <- 1/p_est
    coeffs[ii,] <- ComputeCoeffs(lc,omega_est,tem)
}


### TODO: the des l.c.s do not match fits that well
## plot folded light curve
ii <- 1
par(mfcol=c(2,1))
plotLC(lcs_sdss[[ii]],periods[ii],coeffs[ii,],tem)
plotLC(lcs_des[[ii]],periods[ii],coeffs[ii,],tem)


#####
### add Y band template that is identical to z
### this will be updated in make_template_des.R
# absolute mags
createY <- TRUE
if(createY){
    tem$betas['Y'] <- tem$betas['z']
    tem$betas <- tem$betas[order(names(tem$betas))]
    ## dust
    tem$dust['Y'] <- tem$dust['z']
    tem$dust <- tem$dust[order(names(tem$dust))]
    ## model error
    tem$model_error['Y'] <- tem$model_error['z']
    tem$model_error <- tem$model_error[order(names(tem$model_error))]
    ## templates
    te <- matrix(0,nrow=nrow(tem$templates)+1,ncol=ncol(tem$templates))
    rownames(te) <- c(rownames(tem$templates),"Y")
    te[1:nrow(tem$templates),] <- tem$templates
    te["Y",] <- te["z",]
    te <- te[order(rownames(te)),]
    tem$templates <- te
    ## templatesd
    te <- matrix(0,nrow=nrow(tem$templatesd)+1,ncol=ncol(tem$templatesd))
    rownames(te) <- c(rownames(tem$templatesd),"Y")
    te[1:nrow(tem$templatesd),] <- tem$templatesd
    te["Y",] <- te["z",]
    te <- te[order(rownames(te)),]
    tem$templatesd <- te
    ## tempalate_funcs
    tem$template_funcs["Y"] <- tem$template_funcs["z"]
    tem$template_funcs <- tem$template_funcs[order(names(tem$template_funcs))]
    ## tempalated_funcs
    tem$templated_funcs["Y"] <- tem$templated_funcs["z"]
    tem$templated_funcs <- tem$templated_funcs[order(names(tem$templated_funcs))]
}





## set dust to new values
extc <- read.table("extc.dat",stringsAsFactors=FALSE)
extc <- extc[extc[,1]=="DES",c(2,3)]
tem$dust[extc$V2] <- extc$V3



## how do we empirically estimate new absolute magnitudes?
lcs_resid <- lcs_des
for(ii in 1:length(lcs_des)){
    omega_est <- 1/periods[ii]
    lcs_resid[[ii]][,1] <- (lcs_des[[ii]][,1]*omega_est + coeffs[ii,4]) %% 1.0
    lcs_resid[[ii]][,3] <- lcs_des[[ii]][,3] - PredictTimeBand(lcs_des[[ii]][,1],lcs_des[[ii]][,2],omega_est,coeffs[ii,],tem)
}

lcs_resid <- do.call(rbind,lcs_resid)

bs <- unique(lcs_resid$band)

for(ii in 1:length(bs)){
    lcs_resid_band <- lcs_resid[lcs_resid$band==bs[ii],]
    pdf(paste0("quality_check_des/residuals_",bs[ii],"_new.pdf"))
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
    points(lc_temp[,1],out,type='l',col='red',lwd=2)
    dev.off()
}
