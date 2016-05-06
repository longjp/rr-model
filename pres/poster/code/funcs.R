get_tmss <- function(folder.data){
    folders <- list.dirs(folder.data)
    folders <- folders[folders != folder.data]
    folders <- paste(folders,"/",sep="")
    folders.names <- vapply(strsplit(folders,"/"),
                            function(x){x[length(x)]},c("hello"))
    ## get all light curves in tms form
    fnames <- list()
    for(ii in 1:length(folders)){
        fnames[[ii]] <- list.files(folders[ii])
    }
    fnames <- unique(unlist(fnames))
    tmss <- list()
    for(ii in 1:length(fnames)){
        tmss[[ii]] <- list()
        for(jj in 1:length(folders)){
            fname <- paste(folders[jj],fnames[ii],sep="")
            if(file.exists(fname))  tmss[[ii]][[folders.names[jj]]] <- read.table(fname)
        }
    }
    names(tmss) <- fnames
    return(tmss)
}

plot_lc <- function(lc,period=0,model=rep(0,3),fit=FALSE,fold=FALSE,xlab="Time (Days)"){
    par(mar=c(4.5,4.5,1,.5))
    ylim <- rev(range(lc[,2]))
    shift <- min(lc[,1])
    lc[,1] <- lc[,1] - shift
    if(fold){
        lc[,1] <- lc[,1] %% period
    }
    plot(lc[,1],lc[,2],
         ylim=ylim,
         xlab=xlab,
         ylab="Magnitudes",
         cex.lab=1.5)
    segments(lc[,1],lc[,2] - 2*lc[,3],
             lc[,1],lc[,2] + 2*lc[,3])
    if(fit){
        xs <- seq(min(lc[,1])+shift,max(lc[,1])+shift,length.out=1000)
        ys <- model[1] + model[2]*sin(2*pi*xs/period + model[3])
        points(xs-shift,ys,type='l',col='red')
    }

}

construct_design <- function(w,K,t){
    predesign <- w*outer(t,1:K)
    return(cbind(1,cos(predesign),sin(predesign)))
}

get_freqs <-function(period_min,period_max,freq_del = 0.1/4000){
    freq_max <- 1/period_min
    freq_min <- 1/period_max
    return(2 * pi * seq(freq_min, freq_max, freq_del))
}

compute_params <- function(w,K,mag,weights,X){
    B <- t(X) %*% (X * weights)
    d <- t(X) %*% (mag * weights)
    return(solve(B,d))
}   

compute_rss <- function(w,K,lc){
    X <- construct_design(w,K,lc[,1])
    beta <- compute_params(w,K,lc[,2],lc[,3]^{-2},X)
    r <- (lc[,2] - X%*%beta)
    return(sum(lc[,3] * (r^2)))
}

