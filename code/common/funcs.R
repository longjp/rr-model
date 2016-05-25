within_x <- function(estimates,truth,thres){
    return(apply(abs((estimates - truth)/truth),1,min) < thres)
}

SplitData <- function(tm){
    levs <- levels(tm$band)
    lc <- list()
    for(ii in 1:length(levs)){
        lc[[ii]] <- tm[tm$band==levs[ii],c("time","mag","error")]
    }
    names(lc) <- levs
    return(lc)
}

is_local_min <- function(x){
    x.max <- max(x)
    x.extend <- c(x.max+1,x,x.max+1)
    x.right <- x - x.extend[-c(1,2)]
    x.left <- x - x.extend[1:length(x)]
    return(x.right < 0 & x.left < 0)
}

sort_local_min <- function(x,y){
    loc <- is_local_min(y)
    x <- x[loc]
    y <- y[loc]
    ords <- order(y)
    x <- x[ords]
    return(x)
}

## returns first N elements of x which are
## all separate from each other by t
SeparateBest <- function(x,t,N){
    out <- rep(x[1],N)
    ii <- 2
    while((length(x) >= 2) & (ii <= N)){
        x <- x[2:length(x)]
        x <- x[abs(x - out[ii-1]) > t]
        if(length(x > 0)){
            out[ii] <- x[1]
        }
        ii <- ii + 1
    }
    return(out)
}

## phase shift elements in mu by ix
phase_shift <- function(mu,ix){
    n <- length(mu)
    ix <- ix %% n
    if(ix > 0.5){
        mu <- c(mu[(ix+1):n],mu[1:ix])
    }
    return(mu)
}

phase <- function(x,niter=10){
    mu <- colMeans(x)
    N <- ncol(x)
    n <- nrow(x)
    ## compute denominator, does not depend on del
    xs <- cbind(x[,N],x[,1:(N-1)])
    xd <- x - xs
    den <- rowSums(xd*xd)
    del_total <- rep(0,length(den))
    ## compute numerator, updating mu each time
    for(jj in 1:niter){
        num <- rowSums((t(t(x) - mu))*xd)
        del <- round(-num/den %% N)
        del_total <- (del + del_total) %% N
        x <- t(vapply(1:n,function(ii){phase_shift(x[ii,],del[ii])},rep(0,N)))
        xd <- t(vapply(1:n,function(ii){phase_shift(xd[ii,],del[ii])},rep(0,N)))
        mu <- colMeans(x)
    }
    ## return the total phase shift
    return(del_total)
}
