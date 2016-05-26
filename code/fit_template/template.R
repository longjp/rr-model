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
        if(ix2 >= ix1){
            gammaf[ix1:ix2] <-  temp_funcs[[jj]](t[ix1:ix2])
        }
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

FitTemplate <- function(tm,omegas,tem,NN=1){
    dat <- AugmentData(tm,tem$dust,tem$betas)
    m <- dat[[1]]$mag
    dust <- dat[[1]]$dust
    t <- dat[[1]]$time
    nb <- dat[[2]]
    coeffs <- c(0,0,0,runif(1))
    rss_max <- sum(lm(m~dust)$residuals^2)
    rss <- rep(0,length(omegas))
    for(ii in 1:length(omegas)){
        for(jj in 1:NN){
            coeffs <- NewtonUpdate(coeffs[4],omegas[ii],m,t,dust,nb,tem$template_funcs,tem$templated_funcs)
        }
        gammaf <- ConstructGamma(t,nb,coeffs[4],omegas[ii],tem$template_funcs)
        rss[ii] <- min(sum((m - coeffs[1] - coeffs[2]*dust - coeffs[3]*gammaf)^2),rss_max)
    }
    return(rss)
}

## for a given omega, find coefficients
ComputeCoeffs <- function(tm,omega,tem,NN=10){
    dat <- AugmentData(tm,tem$dust,tem$betas)
    m <- dat[[1]]$mag
    dust <- dat[[1]]$dust
    t <- dat[[1]]$time
    nb <- dat[[2]]
    coeffs <- c(0,0,0,runif(1))
    while(coeffs[3]==0){
        for(jj in 1:NN){
            coeffs <- NewtonUpdate(coeffs[4],omega,m,t,dust,nb,tem$template_funcs,tem$templated_funcs)
        }
    }
    return(coeffs)
}

AmpAlphaDustUpdate <- function(phi,omega,m,t,dust,nb,template_funcs){
    gammaf <- ConstructGamma(t,nb,phi,omega,template_funcs)
    est <- ComputeBeta(m,dust,gammaf)
    alpha <- est["alpha"]
    a <- est["a"]
    d <- est["d"]
    if(a < 0) {
        a <- 0
    }
    out <- c(alpha,d,a)
    names(out) <- NULL
    return(out)
}

ComputeRSSPhase <- function(tm,omega,tem,phis=(1:100)/100){
    dat <- AugmentData(tm,tem$dust,tem$betas)
    m <- dat[[1]]$mag
    dust <- dat[[1]]$dust
    t <- dat[[1]]$time
    nb <- dat[[2]]
    rss_max <- sum(lm(m~dust)$residuals^2)

    rss <- rep(0,length(phis))
    for(ii in 1:length(phis)){
        coeffs <- AmpAlphaDustUpdate(phis[ii],omega,m,t,dust,nb,tem$template_funcs)
        gammaf <- ConstructGamma(t,nb,phis[ii],omega,tem$template_funcs)
        rss[ii] <- min(sum((m - coeffs[1] - coeffs[2]*dust - coeffs[3]*gammaf)^2),rss_max)
    }
    return(rss)
}
