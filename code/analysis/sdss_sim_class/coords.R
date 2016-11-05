rm(list=ls())
library(MASS)
library(rgl)
library(RColorBrewer)
source('coords_funcs.R')
options(width=120)


## need to run this 3x: full data, lomb scargle, and rr lyrae model


###### FOR MAKING "TRUE" MAP
## rr_class <- read.table("../../data/raw/apj326724t2_mrt.txt",skip=42)
## rr_class <- rr_class[,1:2]
## names(rr_class) <- c("ID","class")
## rr <- read.table("../../data/raw/apj326724t3_mrt.txt",skip=30)
## rr <- rr[,1:6]
## names(rr) <- c("ID","ra","dec","ar","d","dg")
## rr <- merge(rr,rr_class)
## rr <- rr[rr$class=="ab",]
## head(rr)
## rr <- rr[,c("ra","dec","d")]


load("rf_rr.RData")
rr <- rf_rr



#### what is input to function
##  rrlyrae: ra, dec, d
##  potting limits


ds <- c(0,120,120,120)
ra <- c(0,20,0,4)*360/24
pts <- EquatToXY(ds,ra)
xlim <- range(pts[,1])
ylim <- range(pts[,2])


rr <- cbind(rr,EquatToXY(rr$d,rr$ra))


plot(rr$x,rr$y,xlim=xlim,ylim=ylim)
##fit.dens <- kde2d(rr$x,rr$y,n=200,lims=c(xlim,ylim),h=2)
fit.dens <- NearestNeighborDensity2d(cbind(rr$x,rr$y),
                                     x1r=xlim,
                                     x2r=ylim,
                                     n=400,k=8)



## convert to cartesian galactocentric coordinates
bls <- EquatorialToGalactic(rr$ra,rr$dec)
cartesian <- SphericalToCartesian(rr$d,bls[,1],bls[,2])
cartesian_mw <- CartesianSunToMW(cartesian)
rr <- cbind(rr,bls)



## make sure this is working right
x <- fit.dens$x
y <- fit.dens$y
dens.mat <- HaloOnCartesian(x,y)



## cdh scales density estimate by volume represented by point
grid.xy <- MakeGrid(x,y)
grid.equat <- XYtoEquat(grid.xy[,1],grid.xy[,2])
grid.ind <- (((grid.equat[,"ra"] <= 360*(4/24)) |
                    (grid.equat[,"ra"] >= 360*(20/24)))
    & (grid.equat[,"d"] <= 120))



cdh <- outer(fit.dens$x,fit.dens$y,
             FUN=function(x,y){return(sqrt(x^2 + y^2))})
grid.ind.mat <- MakeMatrix(grid.ind)
dens.mat[!grid.ind.mat] <- 0
fit.dens$z[!grid.ind.mat] <- 0
dens.mat <- (dens.mat*cdh)  / sum(dens.mat*cdh)
fit.dens$z <- (fit.dens$z*cdh)  / sum(fit.dens$z*cdh)


## does cdh affect halo model too?

z.cont <- fit.dens$z - .8*dens.mat
z.cont <- apply(z.cont,c(1,2),function(x){max(x,0)})


##BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral

MakeContour(z.cont,grid.xy,grid.ind,rr)

#### TODO:
## 1. fix edge effects of density estimator
## 2. clean up code, functionalize
## 3. make plot with downsampled data



