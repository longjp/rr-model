rm(list=ls())
set.seed(1234)

## load necessary libraries
library('parallel')
library('multiband')
library('randomForest')
load("../../fit_template/template_sdss.RData")
source("../../fit_template/fit_template.R")
source("../../common/funcs.R")
source("../funcs.R")


## data source
load("../../data/clean/sdss_sim_class.RData")
load("results.RData")
source("../params.R")


ParameterFits <- function(lcs,period_est){
    N <- length(period_est)
    res <- rep(0,N)
    coeffs <- matrix(0,ncol=4,nrow=N)
    all_res <- vector("list",N)
    for(ii in 1:N){
        omega <- 1/period_est[ii]
        lc <- lcs[[ii]]
        coes <- ComputeCoeffs(lc,omega,tem)
        coeffs[ii,] <- coes
        tm <- LCtoTM(lc)
        dev <- vector("list",length(tm))
        for(jj in 1:length(tm)){
            pred <- (coes[1] + tem$betas[jj] + coes[2]*tem$dust[jj]
                + coes[3]*tem$template_funcs[[jj]]((tm[[jj]][,1]*omega + coes[4]) %% 1))
            dev[[jj]] <- abs((pred - tm[[jj]][,2]))
        }
        all_res[[ii]] <- unlist(dev)
        res[ii] <- median(unlist(dev))
        ##mean_mags <- vapply(tm,function(x){mean(x[,2])},c(0))
    }
    features <- cbind(coeffs,period_est,res)
    colnames(features) <- c("mu","E[B-V]","a","phi","period","res")
    features[,3] <- features[,3]*band_amps['g']
    return(list(features=features,all_res=all_res))
}


period_est_vec <- period_est_FULL[,1]
lcs <- lapply(tms_FULL,function(x){TMtoLC(x)})


band_amps <- (apply(tem$templates,1,max) - apply(tem$templates,1,min))/2
features <- ParameterFits(lcs,period_est_vec)
all_res <- features$all_res
features <- features$features


all_res <- all_res[cl=="rr"]
all_res <- unlist(all_res)
median(all_res)
ct <- quantile(all_res,c(.97))
pdf("absolute_residuals.pdf",width=7,height=6)
hist(all_res[all_res < ct],xlab="Absolute Residual ie |actual mag - model predicted mag|",
     main="RR Lyrae Model for SDSS Stripe 82")
dev.off()
## median residual is 0.03318
##
error <- lapply(lcs,function(x){x$error})
median(abs(rnorm(100000)))*median(unlist(error))
median(all_res)




cols <- c("#00000030",'red')
names(cols) <- c("not","rr")
pdf("residual_versus_amplitude_sdss.pdf")
par(mar=c(5,5,1,1))
plot(features[,'res'],features[,'a'],xlim=c(0,.2),ylim=c(0,1.5),
     col=cols[cl],xlab="median absolute residual",
     ylab="amplitude in g--band",cex.lab=1.5,xaxs='i',yaxs='i')
dev.off()





## same plot for downsampled
period_est_vec <- period_est[,1]
lcs <- lapply(tms,function(x){TMtoLC(x)})


band_amps <- (apply(tem$templates,1,max) - apply(tem$templates,1,min))/2
features <- ParameterFits(lcs,period_est_vec)
all_res <- features$all_res
features <- features$features



cols <- c("#00000030",'red')
names(cols) <- c("not","rr")
pdf("residual_versus_amplitude_sdss_downsampled.pdf")
par(mar=c(5,5,1,1))
plot(features[,'res'],features[,'a'],xlim=c(0,.2),ylim=c(0,1.5),
     col=cols[cl],xlab="median absolute residual",
     ylab="amplitude in g--band",cex.lab=1.5,xaxs='i',yaxs='i')
dev.off()

