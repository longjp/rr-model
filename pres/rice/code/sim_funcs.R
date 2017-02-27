## functions used by sim.R

## compute asymptotic variance, sigma known
Nu1 <- function(A,B,sig,w,n) ((mean(w^2)*A + mean(sig^2*w^2)*B)) / (n*mean(w)^2)

## compute asymptotic variance, sigma known or unknown
Nu2 <- function(B,r,w,x){
    n <- length(x)
    a1 <- mean(r^2*w^2)
    a2 <- mean(r^2*w^2*x)
    a4 <- mean(r^2*w^2*x^2)
    M <- matrix(c(a1,a2,a2,a4),nrow=2)
    return((B%*%M%*%B)/(n*mean(w)^2))
}

## use A,B to compute Delta
ComputeDelta <- function(A,B) max(sum(diag(A))/sum(diag(B)),0)

ComputeAhat <- function(y,x,sigma,Bhat,w=rep(1,length(x)),JJ=2){
    for(jj in 1:JJ){
        lm.fit <- lm(y~x,weights=w)
        gs <- lm.fit$residuals^2 - sigma^2
        gam <- sigma^{-4}
        gamsuminv <- sum(gam)^{-1}
        a1 <- gamsuminv*sum(gs*gam)
        a2 <- gamsuminv*sum(x*gs*gam)
        a4 <- gamsuminv*sum(x^2*gs*gam)
        Ap <- matrix(c(a1,a2,a2,a4),nrow=2)
        Ahat <- Bhat%*%Ap%*%Bhat
        w <- (sigma^2 + ComputeDelta(Ahat,Bhat))^{-1}
    }
    return(Ahat)
}

ComputeBhat <- function(x){
    n <- length(x)
    X <- cbind(1,x)
    return(solve(t(X)%*%X/n))
}

ComputeChat <- function(df){
        a1 <- sum(df[,1])
        a2 <- sum(df[,1]*df[,2])
        a4 <- sum(df[,1]*df[,2]^2)
        return(matrix(c(a1,a2,a2,a4),nrow=2)/nrow(df))
}
    
ComputeUnknownWeights <- function(y,x,m,w=rep(1,length(x)),JJ=2){
    Bhat <- ComputeBhat(x)
    for(jj in 1:JJ){
        res2 <- lm(y~x,weights=w)$residuals^2
        Chats <- by(cbind(res2,x),m,ComputeChat)
        Wk <- vapply(Chats,function(Z){sum(diag(Bhat))/sum(diag(Bhat%*%Z%*%Bhat))},c(0))
        w <- Wk[m]
    }
    return(w)
}
