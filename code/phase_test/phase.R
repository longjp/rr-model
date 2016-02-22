rm(list=ls())
set.seed(1235)

N <- 50
t <- (1:N)/N

sinef <- function(t){
    t <- t %% 1
    return(sin(2*pi*t))
}



n <- 50
x <- matrix(0,nrow=n,ncol=N)
for(ii in 1:n){
    phi <- runif(1,min=0.25,max=.75)
    a <- runif(1,min=0.75,max=1.25)
    x[ii,] <- a*sinef(t + phi)
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
    ## compute numerator, updating mu each time
    for(jj in 1:niter){
        num <- rowSums((t(t(x) - mu))*xd)
        del <- round(-num/den %% N)
        x <- t(vapply(1:n,function(ii){phase_shift(x[ii,],del[ii])},rep(0,N)))
        xd <- t(vapply(1:n,function(ii){phase_shift(xd[ii,],del[ii])},rep(0,N)))
        mu <- colMeans(x)
    }
    return(x)
}


x_new <- phase(x,niter=10)

par(mfcol=c(1,2))
ylim <- range(x)
plot(0,0,col=0,ylim=ylim,xlim=c(0,1),xlab="",ylab="")
for(ii in 1:nrow(x)){
    points(t,x[ii,],type='l',col="#00000050")
}
mu <- colMeans(x)
points(t,mu,type='l',col='red',lwd=2)



ylim <- range(x_new)
plot(0,0,col=0,ylim=ylim,xlim=c(0,1),xlab="",ylab="")
for(ii in 1:nrow(x)){
    points(t,x_new[ii,],type='l',col="#00000050")
}
mu <- colMeans(x_new)
points(t,mu,type='l',col='red',lwd=2)


