rm(list=ls())
library(MASS)
library(rgl)
library(RColorBrewer)
options(width=120)

## TODO: fix XYtoEquat, some issue with arctan

## cartesian_mw has cartesian coordinates
## of RR Lyrae in cartesian coordinates with mw center at origin
XYtoEquat <- function(x,y){
    d <- sqrt(x^2 + y^2)
    ra <- ((atan2(-x,y) + 2*pi) %% (2*pi))*(360/(2*pi))
    return(cbind(d=d,ra=ra,dec=0))
}

EquatToXY <- function(d,ra){
    x <- -d*sin(2*pi*ra/360)
    y <- d*cos(2*pi*ra/360)
    return(cbind(x,y))
}


d <- c(0,120,120,120)
ra <- c(0,20,0,4)*360/24
pts <- EquatToXY(d,ra)
xlim <- range(pts[,1])
ylim <- range(pts[,2])


## ra <- 3*pi/4
## x <- sin(ra)
## y <- cos(ra)
## ra
## x
## y
## atan2(x,y)
## XYtoEquat(x,y)


## rr <- read.table("../../data/raw/apj326724t3_mrt.txt",skip=30)
## rr <- rr[,1:6]
## names(rr) <- c("ID","RA","Decl","ar","dh","dg")
## rr <- cbind(rr,EquatToXY(rr$dh,rr$RA))


## plot(rr$x,rr$y)
## fit.kde <- kde2d(rr$x,rr$y,n=50)



## x <- fit.kde$x
## y <- fit.kde$y
## grid.xy <- cbind(rep(x,length(x)),rep(y,each=length(y)))
## plot(grid.xy,col="#00000030")
## points(rr$x,rr$y,col='red',pch=19)



## grid.equat <- EquatToXY(grid.xy[,1],grid.xy[,2])


## test <- XYtoEquat(rr$x,rr$y)[,1:2]
## head(test)
## plot(test[,1],rr$dh)
## abline(a=0,b=1)

## head(grid.equat)

## grid.equat.ind <- (grid.equat[,"ra"] < 360*(4/24))# | (grid.equat[,"ra"] > 360*(20/24)) 

## par(mfcol=c(1,2))
## plot(-rr$dh*sin(2*pi*rr$RA/360),rr$dh*cos(2*pi*rr$RA/360))
## plot(grid.xy,col=grid.equat.ind)











rr_class <- read.table("../../data/raw/apj326724t2_mrt.txt",skip=42)
rr_class <- rr_class[,1:2]
names(rr_class) <- c("ID","class")



rr <- read.table("../../data/raw/apj326724t3_mrt.txt",skip=30)
rr <- rr[,1:6]
names(rr) <- c("ID","RA","Decl","ar","dh","dg")
rr <- cbind(rr,EquatToXY(rr$dh,rr$RA))
rr <- merge(rr,rr_class)
rr <- rr[rr$class=="ab",]



plot(rr$x,rr$y,xlim=xlim,ylim=ylim)
fit.kde <- kde2d(rr$x,rr$y,n=50,lims=c(xlim,ylim))



plot(rr$x,rr$y,xlab="x",ylab="y",cex.lab=1.3,xlim=xlim,ylim=ylim)
contour(fit.kde,add=TRUE,col="#00000080",cex=2,lims=c(xlim,ylim))



EquatorialToGalactic <- function(alpha,delta){
    alpha <- 2*pi*(alpha/360)
    delta <- 2*pi*(delta/360)
    alpha.ngp <- 2*pi*(192.85 / 360)
    delta.ngp <- 2*pi*((27 + 8/60) / 360)
    alpha.0 <- 2*pi*(282.85 / 360)
    l.0 <- 2*pi*(32.93 / 360)
    b <- (asin(sin(delta)*sin(delta.ngp) -
               cos(delta)*cos(delta.ngp)*sin(alpha-alpha.0)))
    l <- acos(cos(alpha-alpha.0)*cos(delta)/cos(b)) + l.0
    return(cbind(b=(360*b)/(2*pi),l=(360*l)/(2*pi)))
}


SphericalToCartesian <- function(r,b,l){
    phi <- 2*pi*l/360
    theta <- pi/2 - 2*pi*b/360 ## b \in [-90,90] for astronomers, want b \in [0,\pi]
    x <- r*sin(theta)*cos(phi)
    y <- r*sin(theta)*sin(phi)
    z <- r*cos(theta)
    return(cbind(x=x,y=y,z=z))
}

CartesianSunToMW <- function(cartesian){
    SUNtoMWC <- 8 ## distance from sun to milky way center
    cartesian_mw <- t(t(cartesian) - c(SUNtoMWC,0,0))
    return(cartesian_mw)
}


## convert to cartesian galactocentric coordinates
bls <- EquatorialToGalactic(rr$RA,rr$Decl)
cartesian <- SphericalToCartesian(rr$dh,bls[,1],bls[,2])
cartesian_mw <- CartesianSunToMW(cartesian)
rr <- cbind(rr,bls)




## examine results
plot(bls,xlab="b",ylab='l',main="Stripe 82 in Galactic Coordinates")
plot3d(cartesian)
par(mfcol=c(1,3))
plot(rr$dh,sqrt(rowSums(cartesian^2)))
abline(a=0,b=1)
plot(rr$dg,sqrt(rowSums(cartesian^2)))
abline(a=0,b=1)
plot(rr$dg,sqrt(rowSums(cartesian_mw^2)))
abline(a=0,b=1)


MakeGrid <- function(x,y)
    return(cbind(rep(x,length(x)),rep(y,each=length(y))))
}



## make sure this is working right
x <- fit.kde$x
y <- fit.kde$y
grid.xy <- MakeGrid(x,y)
plot(grid.xy,col="#00000030",xlim=xlim,ylim=ylim)
points(rr$x,rr$y,col='red',pch=19)



grid.equat <- XYtoEquat(grid.xy[,1],grid.xy[,2])





head(grid.equat)
grid.equat.ind <- (((grid.equat[,"ra"] < 360*(4/24)) | (grid.equat[,"ra"] > 360*(20/24)))
    & (grid.equat[,"d"] < 120))

par(mfcol=c(1,2))
plot(-rr$dh*sin(2*pi*rr$RA/360),rr$dh*cos(2*pi*rr$RA/360),xlim=xlim,ylim=ylim)
plot(grid.xy,col=grid.equat.ind,xlim=xlim,ylim=ylim)

### but no clear definition for how far we can observe
### so still confused by this halo model calculations


grid.galac <- cbind(grid.equat[,1],EquatorialToGalactic(grid.equat[,2],grid.equat[,3]))
grid.galac_cart <- SphericalToCartesian(grid.galac[,1],grid.galac[,2],grid.galac[,3])
grid.galac_cart_mw <- CartesianSunToMW(grid.galac_cart)
dens <- 1/sqrt(rowSums(grid.galac_cart_mw^2))^3

MakeMatrix <- function(gr){
    return(matrix(gr,nrow=sqrt(length(gr))))
}


fit.kde <- kde2d(rr$x,rr$y,n=50,lims=c(xlim,ylim))
dens.mat <- MakeMatrix(dens)
grid.equat.ind.mat <- MakeMatrix(grid.equat.ind)
dens.mat[!grid.equat.ind.mat] <- 0
fit.kde$z[!grid.equat.ind.mat] <- 0
dens.mat <- dens.mat  / sum(dens.mat)
fit.kde$z <- fit.kde$z  / sum(fit.kde$z)



z.cont <- log10(fit.kde$z / dens.mat)
z.cont[is.na(z.cont)] <- min(z.cont,na.rm=TRUE)



DrawRALine <- function(d){
    ras <- 0:360
    return(EquatToXY(rep(d,length(ras)),ras))
}



##BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral

ncol <- 7
cols <- brewer.pal(ncol,"RdYlBu")
ds <- lapply(25*(1:5),function(x){DrawRALine(x)})
filled.contour(x,y,z.cont,col=cols,nlevels=ncol-3,
               plot.axes = { axis(1);
                   axis(2);
                   points(rr$x,rr$y,col="#00000020");
                   for(ii in ds){points(ii,type='l')}})


plot(rr$x,rr$y)
points(d25,type='l')


head(rr)



### adjust for fact that as distance increases, volume is increasing
### use function better than filled.contour for making contour
### sesar uses subset of these
       ## only rrab stars





## contour(fit.kde,levels=quantile(dens.mat,c(.5,.8,.9)),
##         plot.axes = { axis(1); axis(2); points(x,y,col="#00000020") })




## contour(x, y, dens.mat, col = "pink", method = "edge",
##         vfont = c("sans serif", "plain"),levels=quantile(dens.mat,c(.5,.8,.9)))


## normalize 1/R3 density across area, subtract from observed, plot
