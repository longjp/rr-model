## functions for fitting rr lyrae light curve

## fits RR Lyrae template
##
## arguments
##           lc : light curve, data frame with columns time, band, mag, error
##       omegas : vector of frequencies
##          tem : input templates
##           NN : number of newton steps a each frequency
##   use.errors : should photometric errors be used, generally not advised
##
##
## value
##          rss : the residual sum of squares at each frequency in omegas     
FitTemplate <- function(lc,omegas,tem,NN=5,use.errors=FALSE){
    tem <- CheckTemLC(tem,lc)
    dat <- AugmentData(lc,tem$dust,tem$betas,use.errors)
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

## for a given omega, return cofficients
## example usage: run FitTemplate to determine rss as function of omegas,
## find omega which minimizes rss, then find coeffs for this omega using ComputeCoeffs
##
## arguments
##           lc : light curve, data frame with columns time, band, mag, error
##        omega : frequency
##          tem : input templates
##           NN : number of newton steps, probably 10+ since no warm start
##   use.errors : should photometric errors be used, generally not advised
##
##
## value
##       coeffs : vector of [distance mod,ebv,peak-to-peak g amp,phase]
ComputeCoeffs <- function(lc,omega,tem,NN=10,use.errors=FALSE){
    tem <- CheckTemLC(tem,lc)
    dat <- AugmentData(lc,tem$dust,tem$betas,use.errors)
    m <- dat[[1]]$mag
    dust <- dat[[1]]$dust
    t <- dat[[1]]$time
    nb <- dat[[2]]
    coeffs <- c(0,0,0,runif(1))
    J <- 0
    ## prevent infinite loops with J
    while(coeffs[3]==0 & J < 10){
        for(jj in 1:NN){
            coeffs <- NewtonUpdate(coeffs[4],omega,m,t,dust,nb,tem$template_funcs,tem$templated_funcs)
        }
        J <- J + 1
    }
    ## convert NewtonUpdate "amplitude" to peak-to-peak g band amplitude
    coeffs[3] <- coeffs[3]*diff(range(tem$templates["g",]))    
    return(coeffs)
}


AugmentData <- function(lc,dust,betas,use.errors=FALSE){
    lc <- lc[order(lc$band),]
    nb <- table(lc$band)
    lc$dust <- rep.int(dust,nb)
    lc$mag <- lc$mag - rep.int(betas,nb)
    lc$band <- NULL
    if(!use.errors){
        lc$error <- 1
    }
    return(list(lc=lc,nb=nb))
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

ComputeRSSPhase <- function(lc,omega,tem,phis=(1:100)/100,use.errors=FALSE){
    tem <- CheckTemLC(tem,lc)
    dat <- AugmentData(lc,tem$dust,tem$betas,use.errors)
    m <- dat[[1]]$mag
    dust <- dat[[1]]$dust
    t <- dat[[1]]$time
    nb <- dat[[2]]
    rss_max <- sum(lm(m~dust)$residuals^2)
    rss <- rep(0,length(phis))
    for(ii in 1:length(phis)){
        coeffs <- AmpAlphaDustUpdate(phis[ii],omega,m,t,dust,nb,
                                     tem$template_funcs)
        gammaf <- ConstructGamma(t,nb,phis[ii],omega,tem$template_funcs)
        rss[ii] <- min(sum((m - coeffs[1] - coeffs[2]*dust - coeffs[3]*gammaf)^2),rss_max)
    }
    return(rss)
}

## check and make tem and lc consistent
## if lc has bands not in tem, stop
## if lc has fewer bands then tem, get rid
##    of these bands in tem
CheckTemLC <- function(tem,lc){
    if(prod(unique(lc$band) %in% names(tem$dust)) != 1){
        print("template bands are:")
        print(names(tem$dust))
        print("lc is:")
        print(lc)
        stop("all lc bands must match template names")
    }
    bs <- names(tem$dust)[names(tem$dust) %in% unique(lc$band)]
    if(length(bs) < length(tem$dust)){
        tem$betas <- tem$betas[bs]
        tem$dust <- tem$dust[bs]
        tem$templates <- tem$templates[bs,]
        tem$templatesd <- tem$templatesd[bs,]
        tem$template_funcs <- tem$template_funcs[bs]
        tem$templated_funcs <- tem$templated_funcs[bs]
    }
    return(tem)
}

TBMEtoLC <- function(time,band,mag,error){
    return(data.frame(time,band,mag,error))
}
