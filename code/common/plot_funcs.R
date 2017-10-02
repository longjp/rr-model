plotLC <- function(lc,p_est,coeffs,tem,main=NULL){
    bandpch <- 1:6
    names(bandpch) <- c("u","g","r","i","z","Y")
    bandcol <- c("dodgerblue3","green","red",
                 "mediumorchid1","black","peachpuff4")
    names(bandcol) <- c("u","g","r","i","z","Y")
    lc1 <- lc
    lc1[,1] <- (lc$time %% p_est)/p_est
    lc2 <- lc1
    lc2[,1] <- lc1[,1] + 1
    lc_temp <-rbind(lc1,lc2)
    if(is.null(main)){
        par(mar=c(5,5,1,1))
        main <- ""
    }
    else {
        par(mar=c(5,5,5,1))
    }
    ## obtain magnitude predictions
    ti <- seq(0,p_est,length.out=100)
    ti <- c(ti,ti+p_est)
    m <- PredictAllBand(ti,1/p_est,coeffs,tem)
    ylim <- rev(range(m,lc$mag))
    plot(lc_temp$time,lc_temp$mag,
         col=bandcol[lc_temp$band],pch=bandpch[lc_temp$band],
         ylim=ylim,
         xlab="Phase",ylab="Magnitude",
         xlim=c(0,2),xaxs='i',cex.axis=1.5,cex.lab=1.5,main=main,cex.main=1.5)
    segments(lc_temp$time,
             lc_temp$mag+lc_temp$error,
             lc_temp$time,
             lc_temp$mag-lc_temp$error,col='grey')
    for(ii in 1:ncol(m)){
        points(ti/p_est,m[,ii],type='l',col=bandcol[colnames(m)[ii]])
    }
    ## only make legend for bands used
    bands <- sort(union(colnames(tem$betas),unique(lc$band)))
    bandpch <- bandpch[names(bandpch) %in% bands]
    bandcol <- bandcol[names(bandcol) %in% bands]
    legend("bottomleft",names(bandcol),col=bandcol,pch=bandpch,lty=1)
}
