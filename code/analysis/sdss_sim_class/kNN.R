## kNN density estimator code
MakeGrid <- function(x,y){
    return(cbind(rep(x,length(x)),rep(y,each=length(y))))
}

MakeMatrix <- function(gr){
    return(matrix(gr,nrow=sqrt(length(gr))))
}

NearestNeighborDensity2d <- function(x,x1r,x2r,n=100,k=sqrt(nrow(x))){
    ## construct grid for evaluating density
    x1 <- seq(x1r[1],x1r[2],length.out=n)
    x2 <- seq(x2r[1],x2r[2],length.out=n)
    v <- (x1[2] - x1[1])*(x2[2] - x2[1])
    y <- MakeGrid(x1,x2)
    ## compute and normalize density on grid
    xt <- t(x)
    z <- apply(y,1,function(w){sort(sqrt(colSums((xt-w)^2)))[k]})
    z <- 1/z^2
    z <- (z / sum(z)) / v
    return(list(x=x1,y=x2,z=MakeMatrix(z)))
}

## n <- 10000
## p <- 2

## x <- matrix(rnorm(n*p),ncol=p)
## ## covM <- matrix(c(1,.5,.5,1),nrow=2)
## ## x <- x%*%covM

## d <- NearestNeighborDensity2d(x,x1r=range(x[,1]),x2r=range(x[,2]))



## plot(x,col="#00000020")
## contour(d$x1,d$x2,d$z,add=TRUE,col="red",cex=2)




