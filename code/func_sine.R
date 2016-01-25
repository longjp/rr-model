construct_design <- function(w,K,t){
    predesign <- w*outer(t,1:K)
    return(cbind(1,cos(predesign),sin(predesign)))
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
