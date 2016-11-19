## kNN density estimator code
MakeGrid <- function(x,y){
    return(cbind(rep(x,length(x)),rep(y,each=length(y))))
}

MakeMatrix <- function(gr){
    return(matrix(gr,nrow=sqrt(length(gr))))
}


DrawDCircle <- function(d){
    ras <- 0:360
    return(EquatToXY(rep(d,length(ras)),ras))
}


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

OblateModel <- function(x){
    nh <- 2.77
    qh <- 0.64
    return((1/sqrt(sum(x[1:2]^2 + (x[3]/qh)^2)))^nh)
}

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


## comput halo model rates on a cartesian grid
HaloOnCartesian <- function(x,y){
    grid.xy <- MakeGrid(x,y)
    grid.equat <- XYtoEquat(grid.xy[,1],grid.xy[,2])
    grid.galac <- cbind(grid.equat[,1],
                        EquatorialToGalactic(grid.equat[,2],grid.equat[,3]))
    grid.galac_cart <- SphericalToCartesian(grid.galac[,1],
                                            grid.galac[,2],
                                            grid.galac[,3])
    grid.galac_cart_mw <- CartesianSunToMW(grid.galac_cart)
    dens.oblate <- apply(grid.galac_cart_mw,1,OblateModel)
    return(MakeMatrix(dens.oblate))
}



MakeContour <- function(z.cont,grid.xy,grid.ind,rr,plot_contour=TRUE){
    xlim <- range(grid.xy[grid.ind,1])
    ylim <- range(grid.xy[grid.ind,2])
    xlim <- c(-120,120)
    ylim <- c(0,130)
    plot(0,0,col=0,xaxs='i',yaxs='i',
         xlab="",ylab="",axes=FALSE,cex.lab=4,
         xlim=xlim,ylim=ylim)
    if(plot_contour){
        level <- quantile(log10(z.cont),c(.7,.8,.9,.95,.99))
        ncol <- length(level) + 1
        cols <- rev(brewer.pal(ncol,"RdBu"))
        z.cont <- as.vector(z.cont)
        cols <- rev(brewer.pal(10,name="RdBu"))
        decLocations <- quantile(z.cont[grid.ind],
                                 probs = seq(0.5,0.99,length.out=9),type=4)
        dec <- findInterval(z.cont,c(-Inf,decLocations, Inf))
        points(grid.xy[grid.ind,],col=cols[dec[grid.ind]],pch=20)
    }
    ## locations of RRL
    points(rr$x,rr$y,pch=20)
    ## y=0 line, ie x-axis
    segments(-125,0,125,0,lwd=3)
    ## concentric circles for distance 
    ds <- lapply(25*(1:4),function(x){DrawDCircle(x)})
    for(ii in ds){points(ii,type='l')}
    ## text on x-axis
    at <- c(-100,-75,-50,-25,0,25,50,75,100)
    mtext(text=abs(at),side=1,line=1.7,at=at,cex=4)
    mtext(text="kpc",side=1,line=3.7,at=0,cex=4)
    ## RA lines
    a <- c(270,300,330,0,30,60,90)
    x.from <- rep(0,length(a))
    x.to <- 122*cos(2*pi*(a+90)/(360))
    y.from <- rep(0,length(a))
    y.to <- 122*sin(2*pi*(a+90)/(360))
    segments(x.from,y.from,x.to,y.to)
    ## plot RA labels
    dist <- 126
    locs <- 30*(1:5)
    te <- locs-90
    te[te < 0] <- 360+te[te<0]
    for(ii in 1:length(locs))
        text(dist*cos(locs[ii]*2*pi/360),
             dist*sin(locs[ii]*2*pi/360),
             bquote(.(te[ii])~degree),cex=3)
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
