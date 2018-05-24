## compute uncertainties on all parameter estimates
## output sds (square root diagonal of hessian) along with features
rm(list=ls())


## load data
load("../../data/clean/sdss_sim_class.RData")

ComputeColorMatrix <- function(tms,band_names){
    ## construct feature matrix
    nfeat <- choose(length(band_names),2)
    color_feats <- matrix(0,nrow=length(tms),ncol=2*nfeat)
    color_feats_names <- rep("",ncol(color_feats))
    ## make color names for feature matrix
    cur <- 1
    for(jj in 1:(length(band_names)-1)){
        for(kk in (jj+1):length(band_names)){
            color_feats_names[cur] <- paste0(band_names[jj],band_names[kk])
            color_feats_names[cur+1] <- paste0(band_names[jj],band_names[kk],"_sd")
            cur <- cur + 2
        }
    }
    colnames(color_feats) <- color_feats_names
    ## populate color_feats matrix
    for(ii in 1:length(tms)){
        tm <- tms[[ii]]
        cur <- 1
        for(jj in 1:(length(band_names)-1)){
            for(kk in (jj+1):length(band_names)){
                mag1t <- tm[[band_names[jj]]]$mag
                mag2t <- tm[[band_names[kk]]]$mag
                if(is.null(mag1t) | is.null(mag2t)){
                    color_feats[ii,cur] <- NA
                    color_feats[ii,cur+1] <- NA
                } else { ## will be NA if either mag1t or mag2t has single observation
                    color_feats[ii,cur] <- mean(mag1t) - mean(mag2t)
                    color_feats[ii,cur+1] <- sqrt(var(mag1t)/length(mag1t) + var(mag2t)/length(mag2t))
                }
                cur <- cur + 2
            }
        }
    }
    return(color_feats)
}


band_names <- unique(unlist(lapply(tms,names)))
feats <- ComputeColorMatrix(tms,band_names)
feats_FULL <- ComputeColorMatrix(tms_FULL,band_names)


## some manual checks
head(feats)
head(feats_FULL)
mean(tms[[2]][["g"]]$mag) - mean(tms[[2]][["r"]]$mag)
sqrt(var(tms[[2]][["g"]]$mag)/length(tms[[2]][["g"]]$mag) + var(tms[[2]][["r"]]$mag)/length(tms[[2]][["r"]]$mag))


## get number of observatiosn for poor and well sampled
n <- vapply(tms,function(x){sum(vapply(x,nrow,c(0)))},c(0))
n_FULL <- vapply(tms_FULL,function(x){sum(vapply(x,nrow,c(0)))},c(0))


hist(n)
hist(n_FULL)


save(feats,feats_FULL,n,n_FULL,cl,file="z2-color_feats.RData")
