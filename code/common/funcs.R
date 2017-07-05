within_x <- function(estimates,truth,thres){
    return(apply(abs((estimates - truth)/truth),1,min) < thres)
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


## shifts x to minimize sum((y - mu)^2) where y is shifted x
PhaseGridOne <- function(x,mu){
    n <- length(x)
    out <- vapply(2:n,function(ix){c(x[ix:n],x[1:(ix-1)])},rep(0,n))
    out <- t(cbind(out,x))
    out_diff <- t(t(out) - mu)
    ix <- which.min(rowSums(out_diff^2))
    return(ix)
}
## shifts curves to minimize difference with mu, grid search
## useful to refine estimates from earlier phase alignment algorithms
PhaseGridAll <- function(X){
    mu <- colMeans(X)
    return(apply(X,1,function(x){PhaseGridOne(x,mu)}) %% ncol(X))
}


## solves the optimization problem
##    argmin_{a,Y} \sum_i ||(X[i,,] - a[i]Y||_F^2
##  by iterating updates for a and Y where we constrain ||a||_2 = 1
##  see description.tex for why this is done
##
##
## arguments
##       X : I x T x B array
##       a : vector length I (initial guess)
##       N : number of iterations
##
## value
##     list with elements
##       a : vector length I
##       Y : matrix dim T x B
##
SolveAGamma <- function(X,a=NULL,N=1000){
    I <- dim(X)[1]
    ## initialize a as normalized unit vector
    if(is.null(a)){
        a <- rep(1,I)
    }
    a <- a / sqrt(sum(a*a))
    ## update a and Y iteratively, N times
    for(jj in 1:N){
        Y <- vapply(1:I,function(ii){a[ii]*X[ii,,]},matrix(0,nrow=dim(X)[2],ncol=dim(X)[3]))
        Y <- apply(Y,1:2,sum)
        a <- vapply(1:I,function(ii){sum(Y*X[ii,,])},c(0))
        a <- a / sqrt(sum(a*a))
    }
    return(list(a=a,Y=Y))
}


LCtoTM <- function(lc){
    lc[,1] <- lc[,1] - min(lc[,1])
    levs <- unique(lc$band)
    levs <- levs[order(levs)]
    tm <- list()
    for(ii in 1:length(levs)){
        tm[[ii]] <- lc[lc$band == levs[ii],c("time","mag","error")]
        names(tm[[ii]]) <- c("time","mag","error")
    }
    names(tm) <- levs
    nb <- vapply(tm,nrow,c(0))
    tm <- tm[nb>0]
    return(tm)
}

TMtoLC <- function(tm){
    nb <- vapply(tm,nrow,c(0))
    bs <- rep.int(names(tm),nb)
    lc <- do.call(rbind,tm)
    lc <- data.frame(time=lc[,1],band=bs,mag=lc[,2],error=lc[,3])
    lc[,2] <- as.character(lc[,2])
    return(lc)
}

    
