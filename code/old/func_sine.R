construct_design <- function(w,K,t){
    predesign <- w*outer(t,1:K)
    return(cbind(1,cos(predesign),sin(predesign)))
}

compute_params <- function(w,K,mag,X){
    B <- t(X) %*% X
    d <- t(X) %*% mag
    return(solve(B,d))
}   

compute_rss <- function(w,K,lc){
    X <- construct_design(w,K,lc[,1])
    beta <- compute_params(w,K,lc[,2],X)
    return(sum((lc[,2] - X%*%beta)^2))
}


get_sinusoidal_params <- function(beta){
    beta0 <- beta[1]
    amp <- rep(0,(length(beta)-1)/2)
    rho <- rep(0,(length(beta)-1)/2)
    for(ii in 1:length(amp)){
        amp[ii] <- sqrt(beta[2*ii]^2 + beta[2*ii + 1]^2)
        rho[ii] <- atan2(beta[2*ii],beta[2*ii+1])
    }
    return(list(beta0=beta0,amp=amp,rho=rho))
}
