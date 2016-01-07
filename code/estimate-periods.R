rm(list=ls())

load("tms.RData")

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


p_est <- rep(0,length(tms))
K <- 1 ## try other K as well
period_min <- 0.2 ## since this is a cepheid, should have period in [1,100]
period_max <- 1
for(ii in 1:length(tms)){
    print(ii)
    lc <- tms[[ii]]
    omegas <- get_freqs(period_min,period_max,.1/diff(range(lc[,1])))
    rss <- vapply(omegas,compute_rss,c(0),K,lc)
    p_est[ii] <- (2*pi)/omegas[which.min(rss)]
}



save(tms,periods,p_est,file="estimates.RData")
