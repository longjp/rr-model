rm(list=ls())
library(MASS)
library(rgl)
library(RColorBrewer)
source('kNN.R')
options(width=120)



rr_class <- read.table("../../data/raw/apj326724t2_mrt.txt",skip=42)
rr_class <- rr_class[,1:2]
names(rr_class) <- c("ID","class")



rr <- read.table("../../data/raw/apj326724t3_mrt.txt",skip=30)
rr <- rr[,1:6]
names(rr) <- c("ID","RA","Decl","ar","dh","dg")
rr <- merge(rr,rr_class)
rr <- rr[rr$class=="ab",]

head(rr)
rr <- rr[,c("RA","Decl","dh")]




#### what is input to function
##  rrlyrae: ra, dec, d
##  potting limits


d <- c(0,120,120,120)
ra <- c(0,20,0,4)*360/24
pts <- EquatToXY(d,ra)
xlim <- range(pts[,1])
ylim <- range(pts[,2])



rr <- cbind(rr,EquatToXY(rr$dh,rr$RA))


plot(rr$x,rr$y,xlim=xlim,ylim=ylim)
##fit.dens <- kde2d(rr$x,rr$y,n=200,lims=c(xlim,ylim),h=2)
fit.dens <- NearestNeighborDensity2d(cbind(rr$x,rr$y),
                                     x1r=xlim,
                                     x2r=ylim,
                                     n=400,k=8)



## convert to cartesian galactocentric coordinates
bls <- EquatorialToGalactic(rr$RA,rr$Decl)
cartesian <- SphericalToCartesian(rr$dh,bls[,1],bls[,2])
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

a <- c(270,300,330,0,30,60,90)
x.from <- rep(0,length(a))
x.to <- 200*cos(2*pi*(a+90)/(360))
y.from <- rep(0,length(a))
y.to <- 200*sin(2*pi*(a+90)/(360))

##BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral

## z.cont <- log10(z.cont)
## z.cont[is.na(z.cont)] <- min(z.cont,na.rm=TRUE)







level <- quantile(log10(z.cont),c(.7,.8,.9,.95,.99))
ncol <- length(level) + 1
cols <- rev(brewer.pal(ncol,"RdBu"))


z.cont <- as.vector(z.cont)



cols <- rev(brewer.pal(10,name="RdBu"))
decLocations <- quantile(z.cont[grid.ind],
                         probs = seq(0.5,0.99,length.out=9),type=4)
dec <- findInterval(z.cont,c(-Inf,decLocations, Inf))


plot(grid.xy[grid.ind,],col=cols[dec[grid.ind]],pch=20,xaxs='i',yaxs='i')
ds <- lapply(25*(1:5),function(x){DrawDCircle(x)})
points(rr$x,rr$y,pch=20)
for(ii in ds){points(ii,type='l')}
segments(x.from,y.from,x.to,y.to)


#### TODO:
## 1. fix edge effects of density estimator
## 2. clean up code, functionalize
## 3. make plot with downsampled data



