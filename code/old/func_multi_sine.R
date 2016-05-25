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


to_list_lc <- function(lc){
    lc.list <- list()
    bands <- unique(lc$ampj)
    for(bb in 1:length(bands)){
        lc.list[[bb]] <- lc[lc$ampj==bands[bb],1:3]
    }
    return(lc.list)
}
