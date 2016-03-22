## functions for fitting rr lyrae light curve

AugmentData <- function(tm,dust,betas,use.errors=FALSE){
    tm <- tm[order(tm$band),]
    nb <- table(tm$band)
    tm$dust <- rep.int(dust,nb)
    tm$mag <- tm$mag - rep.int(betas,nb)
    tm$band <- NULL
    if(!use.errors){
        tm$error <- 1
    }
    return(list(tm=tm,nb=nb))
}

ConstructGamma <- function(t,nb,phi,omega,temp_funcs){
    t <- (t*omega + phi) %% 1
    bix <- c(0,cumsum(nb))
    gammaf <- rep(0,length(t))
    for(jj in 1:(length(bix)-1)){
        ix1 <- (bix[jj]+1)
        ix2 <- bix[jj+1]
        gammaf[ix1:ix2] <-  temp_funcs[[jj]](t[ix1:ix2])
    }
    return(gammaf)
}

## computes beta
ComputeBeta <- function(m,dust,gammaf){
    X <- cbind(alpha=1,d=dust,a=gammaf)
    B <- t(X)%*%X
    d <- t(X)%*%m
    z <- solve(B,d)
   return(z[,1])
}

NewtonUpdate <- function(phi,omega,m,t,dust,nb,template_funcs,templated_funcs){
    gammaf <- ConstructGamma(t,nb,phi,omega,template_funcs)
    est <- ComputeBeta(m,dust,gammaf)
    alpha <- est["alpha"]
    a <- est["a"]
    d <- est["d"]
    if(a > 0){
        gammafd <- ConstructGamma(t,nb,phi,omega,templated_funcs)
        mp <- m - alpha - d*dust
        del <- sum(gammafd*(mp-a*gammaf))
        h <- a*sum(gammafd*gammafd)
        phi <- (phi + h^{-1}*del) %% 1
    } else {
        a <- 0
        phi <- runif(1)
    }
    out <- c(alpha,d,a,phi=phi)
    names(out) <- NULL
    return(out)
}
