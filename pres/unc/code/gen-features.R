rm(list=ls())
library(moments)
library(multiband)
library(parallel)

MC.CORES <- 7

folders <- list.files('data/')
lcs <- vector("list",length(folders))
names(lcs) <- folders
for(ii in 1:length(folders)){
    fnames <- list.files(paste("data/",folders[ii],sep=""),full.names=TRUE)
    lcs[[ii]] <- vector("list",length(fnames))
    names(lcs[[ii]]) <- fnames
    for(jj in 1:length(lcs[[ii]])){
        lcs[[ii]][[jj]] <- as.matrix(read.table(fnames[jj],header=FALSE))
    }
}



### feature functions
quantile95_5 <- function(lc){
    high <- quantile(lc[,2],.95)
    low <- quantile(lc[,2],.05)
    return(as.vector(high - low))
}
skew <- function(lc){
    return(skewness(lc[,2]))
}
period <- function(lc){
    out <- pgls(list(lc),period_min=.2,period_max=500,BCD_flag=FALSE)
    return(out$period_seq_all[which.min(out$rss_ls)])
}



compute_features <- function(lc){
    s <- skew(lc)
    amp <- quantile95_5(lc)
    p <- period(lc)
    ## p2p scatter
    ords <- order(lc[,1])
    lc <- lc[ords,]
    med.all <- mad(lc[,2])
    med.adj <- median(abs(lc[2:nrow(lc),2] - lc[1:(nrow(lc)-1),2]))
    p2p_scatter <- med.adj / med.all
    return(c(p=p,s=s,amp=amp,p2p_scatter=p2p_scatter))
}

features <- vector("list",length(lcs))
names(features) <- names(lcs)
for(ii in 1:length(features)){
    print(paste("computing",names(features)[[ii]]))
    features[[ii]] <- do.call(rbind,mclapply(lcs[[ii]],compute_features,mc.cores=MC.CORES))
    rownames(features[[ii]]) <- names(lcs[[ii]])
}

## merge feature files together to create single output
features <- do.call(rbind,features)
features <- as.data.frame(features)
num_lcs <- vapply(lcs,length,c(0))
features$class <- rep(names(num_lcs),times=num_lcs)
save(features,file="features.RData")
