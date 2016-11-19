rm(list=ls())
library(MASS)
library(RColorBrewer)
source('map_funcs.R')
options(width=120)

########## TO STORE TRUE, RR MODEL, AND LOMB DATA
rr <- vector("list",2)
names(rr) <- c("true","rr_model")

#### FOR MAKING "TRUE" MAP
rr_class <- read.table("../../data/raw/apj326724t2_mrt.txt",skip=42)
rr_class <- rr_class[,1:2]
names(rr_class) <- c("ID","class")
tr <- read.table("../../data/raw/apj326724t3_mrt.txt",skip=30)
tr <- tr[,1:6]
names(tr) <- c("ID","ra","dec","ar","d","dg")
tr <- merge(tr,rr_class)
tr <- tr[tr$class=="ab",]
tr <- tr[,c("ra","dec","d")]
rr[[1]] <- tr

#### FOR MAKING RR MODEL ESTIMATED MAP
load("rf_rr.RData")
rr[[2]] <- rf_rr

#### FOR MAKING LOMB BASED MAP
##  todo


PrepPlot <- function(rr){
    ds <- c(0,120,120,120)
    ra <- c(0,20,0,4)*360/24
    pts <- EquatToXY(ds,ra)
    xlim <- range(pts[,1])
    ylim <- range(pts[,2])

    xlim <- c(-125,125)
    ylim <- c(0,125)
    
    rr <- cbind(rr,EquatToXY(rr$d,rr$ra))
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
                  (grid.equat[,"ra"] >= 360*(20/24 + 32/(24*60))))
        & (grid.equat[,"d"] <= 120)
        & (grid.equat[,"d"] >= 5))
    cdh <- outer(fit.dens$x,fit.dens$y,
                 FUN=function(x,y){return(sqrt(x^2 + y^2))})
    grid.ind.mat <- MakeMatrix(grid.ind)
    dens.mat[!grid.ind.mat] <- 0
    fit.dens$z[!grid.ind.mat] <- 0
    dens.mat <- (dens.mat*cdh)  / sum(dens.mat*cdh)
    fit.dens$z <- (fit.dens$z*cdh)  / sum(fit.dens$z*cdh)
    ## subtract off halo model
    z.cont <- fit.dens$z - .8*dens.mat
    z.cont <- apply(z.cont,c(1,2),function(x){max(x,0)})
    return(list(z.cont=z.cont,grid.xy=grid.xy,grid.ind=grid.ind,rr=rr))

}



## plot with "truth"
out <- PrepPlot(rr[[1]])
png("density_true.png",width = 1400, height = 800, units = "px")
MakeContour(out$z.cont,out$grid.xy,out$grid.ind,out$rr)
dev.off()
png("density_true_points.png",width = 1400, height = 800, units = "px")
MakeContour(out$z.cont,out$grid.xy,out$grid.ind,out$rr,plot_contour=FALSE)
dev.off()

## plot with "estimates"
out <- PrepPlot(rr[[2]])
png("density_rr_model_sampled.png",width = 1400, height = 800, units = "px")
MakeContour(out$z.cont,out$grid.xy,out$grid.ind,out$rr)
dev.off()
png("density_rr_model_sampled_points.png",width = 1400, height = 800, units = "px")
MakeContour(out$z.cont,out$grid.xy,out$grid.ind,out$rr,plot_contour=FALSE)
dev.off()


source('map_funcs.R')
MakeContour(out$z.cont,out$grid.xy,out$grid.ind,out$rr)


#### TODO:
## 1. get rid of RR past 120 kpc and less than 5 kpc
## 2. does cdh affect halo model too?
## 3. fix edge effects of density estimator


