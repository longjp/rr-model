rm(list=ls())
library(MASS)
library(rgl)
options(width=120)

rr <- read.table("../../data/raw/apj326724t3_mrt.txt",skip=30)
names(rr)[1:6] <- c("ID","RA","Decl","ar","dh","dg")

plot(-rr$dh*sin(2*pi*rr$RA/360),rr$dh*cos(2*pi*rr$RA/360))

fit.kde <- kde2d(-rr$dh*sin(2*pi*rr$RA/360),rr$dh*cos(2*pi*rr$RA/360),n=50)


fit.kde$z[20:30,20:30] <- max(fit.kde$z)

plot(-rr$dh*sin(2*pi*rr$RA/360),rr$dh*cos(2*pi*rr$RA/360),
     xlab="x",ylab="y",cex.lab=1.3)
contour(fit.kde,add=TRUE,col="#00000020",cex=2)



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








## cartesian_mw has cartesian coordinates
## of RR Lyrae in cartesian coordinates with mw center at origin
ToPolar <- function(x,y){
    d <- sqrt(x^2 + y^2)
    ra <- (atan2(-x,y)/(2*pi))*360
    return(cbind(d=d,ra=ra,dec=0))
}




## make sure this is working right
x <- fit.kde$x
y <- fit.kde$y
grid.xy <- cbind(rep(x,length(x)),rep(y,each=length(y)))
grid.equat <- ToPolar(grid.xy[,1],grid.xy[,2])


head(grid.equat)
grid.equat.ind <- (grid.equat[,"ra"] < 360*(4/24))# | (grid.equat[,"ra"] > 360*(20/24)) 

par(mfcol=c(1,2))
plot(-rr$dh*sin(2*pi*rr$RA/360),rr$dh*cos(2*pi*rr$RA/360))
plot(grid.xy,col=grid.equat.ind)


grid.galac <- cbind(grid.equat[,1],EquatorialToGalactic(grid.equat[,2],grid.equat[,3]))
grid.galac_cart <- SphericalToCartesian(grid.galac[,1],grid.galac[,2],grid.galac[,3])
grid.galac_cart_mw <- CartesianSunToMW(grid.galac_cart)
dens <- 1/sqrt(rowSums(grid.galac_cart_mw^2))^3
matrix(dens

## normalize 1/R3 density across area, subtract from observed, plot
