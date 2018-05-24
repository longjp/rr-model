### TODO: uncertainties should decrease with increasing n
###       this is not the case for RSS uncertainty


####
#### - There are strong correlations between class and feature errors
#### - larger feature errors suggest not rrlyrae

## compute uncertainties on all parameter estimates
## output sds (square root diagonal of hessian) along with features
rm(list=ls())

require(numDeriv)
load("../../fit_template/template_sdss.RData")
load("../../fit_template/template_des.RData")
source("../../fit_template/fit_template.R")
load("../../fit_template/template_sdss.RData")
source("../funcs.R")
source("../../common/funcs.R")

## load data
load("../../data/clean/sdss_sim_class.RData")
load("z0-fit.RData")


## computes the fisher information, the hessian of the neg log like
## can invert resulting matrix to determine asymptotic variance
##
## arguments
##      coeffs : [mu,ebv,a,phase]
##       omega : frequency
##          lc : light curve
ComputeFI <- function(coeffs,omega,lc,tem,use.errors=TRUE,use.dust=TRUE){
    if(use.dust){
        use.dust <- CheckNumberBands(lc)
    }
    CheckLC(lc)
    ## unpack coefficients
    tem <- CheckTemLC(tem,lc)
    dat <- AugmentData(lc,tem,use.errors)
    mag <- dat[[1]]$mag
    dust <- dat[[1]]$dust
    t <- dat[[1]]$time
    weights <- 1 / dat[[1]]$error^2
    nb <- dat[[2]]
    coeffs <- c(omega,coeffs)
    return(hessian(ComputeRSS,coeffs,mag=mag,tem=tem,nb=nb,t=t,dust=dust,weights=weights))
}

ComputeRSS <- function(coeffs,mag,tem,nb,t,dust,weights){
    omega <- coeffs[1]
    phi <- coeffs[5]
    coeffs <- coeffs[2:4]
    mag <- mag - rep.int(tem$abs_mag(1/omega,tem)[1,],nb)
    gammaf <- ConstructGamma(t,nb,phi,omega,tem$template_funcs)
    resid <- mag - coeffs[1] - coeffs[2]*dust - coeffs[3]*gammaf
    return(sum(weights*resid^2))
}



####### CHECK PERIODS REASONABLE
N <- sum(cl=="rr")
## fraction of times period is best
print("accuracies, top period:")
period_est <- period_est[,1] ## just use best fit period
print(paste("1%:",mean(abs((periods[1:N] - period_est[1:N])/periods[1:N]) < 0.01)))
print(paste("0.1%:",mean(abs((periods[1:N] - period_est[1:N])/periods[1:N]) < 0.001)))
print(paste("0.01%:",mean(abs((periods[1:N] - period_est[1:N])/periods[1:N]) < 0.0001)))

## fraction of times period is best on FULL
print("accuracies, top period:")
period_est_FULL <- period_est_FULL[,1] ## just use best fit period
print(paste("1%:",mean(abs((periods[1:N] - period_est_FULL[1:N])/periods[1:N]) < 0.01)))
print(paste("0.1%:",mean(abs((periods[1:N] - period_est_FULL[1:N])/periods[1:N]) < 0.001)))
print(paste("0.01%:",mean(abs((periods[1:N] - period_est_FULL[1:N])/periods[1:N]) < 0.0001)))




####### ESTIMATE MODEL PARAMETERS AND UNCERTAINTIES
## a,ebv,phase,period,rss,variability measures,colors
phis <- (1:1000)/1000

## on full light curves, get features and uncertainties
coeffs_FULL <- matrix(0,nrow=length(tms_FULL),ncol=4)
for(ii in 1:length(tms_FULL)){
    tm <- tms_FULL[[ii]]
    lc <- TMtoLC(tm)
    p_est <- period_est_FULL[ii]
    omega <- 1 / p_est
    rssphi <- ComputeRSSPhase(lc,omega,tem_sdss,phis)
    phi <- phis[which.min(rssphi)]
    coeffs_FULL[ii,] <- ComputeCoeffsPhase(lc,omega,phi,tem_sdss)
}

coeffsu_FULL <- matrix(0,nrow=length(tms_FULL),ncol=4)
for(ii in 1:nrow(coeffsu_FULL)){
    fi <- ComputeFI(coeffs_FULL[ii,],1/period_est_FULL[ii],TMtoLC(tms_FULL[[ii]]),tem_sdss)
    coeffsu_FULL[ii,] <- diag(solve(fi[2:5,2:5]))
}

## on poorly observed light curves, get features and uncertainties
coeffs <- matrix(0,nrow=length(tms),ncol=4)
for(ii in 1:length(tms)){
    tm <- tms[[ii]]
    lc <- TMtoLC(tm)
    p_est <- period_est[ii]
    omega <- 1 / p_est
    rssphi <- ComputeRSSPhase(lc,omega,tem_sdss,phis)
    phi <- phis[which.min(rssphi)]
    coeffs[ii,] <- ComputeCoeffsPhase(lc,omega,phi,tem_sdss)
}

coeffsu <- matrix(0,nrow=length(tms),ncol=4)
for(ii in 1:nrow(coeffsu)){
    fi <- ComputeFI(coeffs[ii,],1/period_est[ii],TMtoLC(tms[[ii]]),tem_sdss)
    coeffsu[ii,] <- diag(solve(fi[2:5,2:5]))
}


## bind period estimates and coefficients
coeffs <- cbind(period_est,coeffs)
coeffs_FULL <- cbind(period_est_FULL,coeffs_FULL)

## take square root of coeffsu
coeffsu <- sqrt(coeffsu)
coeffsu_FULL <- sqrt(coeffsu_FULL)

## add names to uncertainties
fnames <- c("period","mu","ebv","amp","phase")
fnames_sig <- c("mu_sig","ebv_sig","amp_sig","phase_sig")
colnames(coeffs) <- fnames
colnames(coeffs_FULL) <- fnames
colnames(coeffsu) <- fnames_sig
colnames(coeffsu_FULL) <- fnames_sig






####
#### compute rss and standard deviation of rss
####

ComputeResiduals <- function(lc,omega,coeffs,tem,use.errors,use.dust){
    if(use.dust){
        use.dust <- CheckNumberBands(lc)
    }
    CheckLC(lc)
    tem <- CheckTemLC(tem,lc)
    dat <- AugmentData(lc,tem,use.errors)
    m <- dat[[1]]$mag
    dust <- dat[[1]]$dust
    t <- dat[[1]]$time
    weights <- 1 / dat[[1]]$error^2
    nb <- dat[[2]]
    m <- m - rep.int(tem$abs_mag(1/omega,tem)[1,],nb)
    rss_max <- sum((lm(m~dust,weights=weights)$residuals^2)*weights)
    rss <- rep(0,length(phis))
    gammaf <- ConstructGamma(t,nb,coeffs[4],omega,tem$template_funcs)
    resid <- m - coeffs[1] - coeffs[2]*dust - coeffs[3]*gammaf
    return(weights*resid^2)
}

mss <- matrix(0,nrow=length(tms),ncol=2)
colnames(mss) <- c("mss","mss_sd")
for(ii in 1:length(tms)){
    res <- ComputeResiduals(TMtoLC(tms[[ii]]),1/coeffs[ii,1],coeffs[ii,2:5],
                            tem_sdss,use.errors=TRUE,use.dust=TRUE)
    mss[ii,1] <- sum(res) / (length(res)-9)
    mss[ii,2] <- sd(res) / sqrt(length(res))
}

mss_FULL <- matrix(0,nrow=length(tms_FULL),ncol=2)
colnames(mss_FULL) <- c("mss","mss_sd")
for(ii in 1:length(tms_FULL)){
    res <- ComputeResiduals(TMtoLC(tms_FULL[[ii]]),1/coeffs_FULL[ii,1],coeffs_FULL[ii,2:5],
                            tem_sdss,use.errors=TRUE,use.dust=TRUE)
    mss_FULL[ii,1] <- sum(res) / (length(res)-9)
    mss_FULL[ii,2] <- sd(res) / sqrt(length(res))
}


npoints <- vapply(lapply(tms,TMtoLC),nrow,c(0))
out <- aggregate(mss[,1],by=list(npoints),FUN=median)
out2 <- aggregate(mss[,2],by=list(npoints),FUN=median)
par(mfcol=c(2,1))
plot(out[,1],out[,2])
plot(out2[,1],out2[,2])


dat <- data.frame(cl=cl,mss)
g <- ggplot(dat, aes(mss))
g + geom_density(aes(fill=factor(cl)), alpha=0.8) + 
    labs(title="Density plot", 
         x="mss",
         fill="class") + xlim(0, 10)

dev.new()
dat <- data.frame(cl=cl,mss_FULL)
g <- ggplot(dat, aes(mss))
g + geom_density(aes(fill=factor(cl)), alpha=0.8) + 
    labs(title="Density plot", 
         x="mss",
         fill="class") + xlim(0, 10)

dat <- data.frame(cl=cl,mss,coeffs)
dat$g <- dat$amp / dat$mss
g <- ggplot(dat, aes(g))
g + geom_density(aes(fill=factor(cl)), alpha=0.8) + 
    labs(title="Density plot", 
         x="g",
         fill="class") + xlim(0, 2)


dev.new()
dat <- data.frame(cl=cl,mss_FULL,coeffs_FULL)
dat$g <- dat$amp / dat$mss
g <- ggplot(dat, aes(g))
g + geom_density(aes(fill=factor(cl)), alpha=0.8) + 
    labs(title="Density plot", 
         x="g",
         fill="class") + xlim(0, 2)
#### TODO: QUESTION: WHY is sd MSS LOWER FOR poorly sampled than well sampled
#### TODO: MSS lower for FULL 
mean(mss[,1] > mss_FULL[,1])
sum(mss_FULL[,2] < mss[,2])












####
#### compute period uncertainty due to vuong
####

## TODO: need to keep 2nd best period estimate for this part

## estimates confidence of period estimate using baluev's method, one-sided
PeriodConfidence <- function(lc,omega1,omega2,tem,phis=(1:100)/100,use.errors=TRUE,use.dust=TRUE){
    phi1 <- phis[which.min(ComputeRSSPhase(lc,omega1,tem,phis,use.errors,use.dust))]
    phi2 <- phis[which.min(ComputeRSSPhase(lc,omega2,tem,phis,use.errors,use.dust))]
    coeffs1 <- ComputeCoeffsPhase(lc,omega1,phi1,tem,use.errors,use.dust)
    coeffs2 <- ComputeCoeffsPhase(lc,omega2,phi2,tem,use.errors,use.dust)
    r1 <- ComputeResiduals(lc,omega1,coeffs1,tem,use.errors,use.dust)
    r2 <- ComputeResiduals(lc,omega2,coeffs2,tem,use.errors,use.dust)
    return(pnorm(0,mean=mean(r1-r2),sd=sd(r1-r2)/sqrt(length(r1))))
}



## TODO: check that probability correct increases with increasing n


####
#### compute colors and uncertainties
####
ComputeColors <- function(tms,color_names){
    com <- t(combn(1:length(color_names), 2))
    cnames <- vapply(1:nrow(com),function(ix){paste(color_names[com[ix,]],collapse="-")},c("0"))
    ## create matrices for storing colors and uncertainties
    cols <- matrix(0,nrow=length(tms),ncol=nrow(com))
    colsu <- matrix(0,nrow=length(tms),ncol=nrow(com))
    colnames(cols) <- cnames
    colnames(colsu) <- cnames
    ## run through all light curves and all filters
    for(ii in 1:length(tms)){
        for(jj in 1:nrow(com)){
            t1 <- tms[[ii]][[color_names[com[jj,1]]]][,2]
            t2 <- tms[[ii]][[color_names[com[jj,2]]]][,2]
            ## NAs if either filter has 0 obs
            if(is.null(t1) | is.null(t2)){
                cols[ii,jj] <- NA
                colsu[ii,jj] <- NA
            } else if((length(t1)==1) | (length(t2)==1)){
                ## NAs if either filter has 1 obs
                cols[ii,jj] <- mean(t1) - mean(t2)
                colsu[ii,jj] <- NA
            } else {
                ## if at least two obs in each filter, compute means and sd
                cols[ii,jj] <- mean(t1) - mean(t2)
                colsu[ii,jj] <- sqrt(var(t1) / (length(t1) - 1)  + var(t2) / (length(t2) - 1))
            }
        }
    }
    return(list(cols=cols,colsu=colsu))
}



### TODO: how well does npoints predict uncertainty
color_names <- unique(unlist(lapply(tms,names)))
temp1 <- ComputeColors(tms,color_names)
temp2 <- ComputeColors(tms_FULL,color_names)






## plot noisy versus clean features
plot(coeffs[,1],coeffs[,4],ylim=c(0,3))
segments(coeffs[,1] - coeffsu[,1],coeffs[,4],coeffs[,1] + coeffsu[,1],coeffs[,4],col="#00000020")
segments(coeffs[,1],coeffs[,4] - coeffsu[,4],coeffs[,1],coeffs[,4] + coeffsu[,4],col="#00000020")
dev.new()
plot(coeffs_FULL[,1],coeffs_FULL[,4],ylim=c(0,3))
segments(coeffs_FULL[,1] - coeffsu_FULL[,1],coeffs_FULL[,4],coeffs_FULL[,1] + coeffsu_FULL[,1],coeffs_FULL[,4],col="#00000020")
segments(coeffs_FULL[,1],coeffs_FULL[,4] - coeffsu_FULL[,4],coeffs_FULL[,1],coeffs_FULL[,4] + coeffsu_FULL[,4],col="#00000020")



temp <- rep("",length(npoints))
temp[npoints < median(npoints)] <- "poor"
clpoor <- paste0(cl,temp)



cols <- c("black","red","orange","blue")
names(cols) <- c("rr","not","rrpoor","notpoor")

pchs <- c(1:4)
names(pchs) <- c("rr","not","rrpoor","notpoor")



## plot noisy versus clean features
plot(coeffs[,1],coeffs[,4],ylim=c(0,3),col=cols[clpoor],pch=pchs[cl])
dev.new()
plot(coeffs_FULL[,1],coeffs_FULL[,4],ylim=c(0,3),col=cols[cl],pch=pchs[cl])


to_use <- npoints >= median(npoints)
plot(coeffs[to_use,1],coeffs[to_use,4],ylim=c(0,3),col=cols[cl[to_use]],pch=pchs[cl[to_use]])
dev.new()
to_use <- npoints >= median(npoints)
plot(coeffs[!to_use,1],coeffs[!to_use,4],ylim=c(0,3),col=cols[cl[!to_use]],pch=pchs[cl[!to_use]])






## 





## read in des data for same light curve
fname <- "LC_402316_des.dat"
lc_des <- read.table(fname,header=TRUE,stringsAsFactors=FALSE)
lc_des <- lc_des[,c(1,4,2,3)]
names(lc_des) <- c("time","band","mag","error")

## construct frequency grid
freq_min <- 1.0/0.95 ## max period of rrab is about 0.95 days 
freq_max <- 1.0/0.4 ## min period of rrab is about 0.4 days
max_phase_error <- .01 ## maximum phase error fraction, see documentation
freq_space <- (max_phase_error*4)/(max(lc_des[,1]) - min(lc_des[,1]))
omegas <- seq(from=freq_min,to=freq_max,by=freq_space)


## estimate period with des templates
rss_des <- FitTemplate(lc_des,omegas,tem_des)
omega_est_des <- omegas[which.min(rss_des)]
p_est_des <- 1/omega_est_des
coeffs_des <- ComputeCoeffs(lc_des,omega_est_des,tem_des)
names(coeffs_des) <- c("mu","ebv","amp","phase")


des_poor <- ComputeFI(coeffs_des,omega_est_des,lc_des,tem_des)
des_poor
