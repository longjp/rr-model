plotLC <- function(lc,p_est,coeffs,tem){
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
    par(mar=c(5,5,1,1))
    plot(lc_temp$time,lc_temp$mag,
         col=bandcol[lc_temp$band],pch=bandpch[lc_temp$band],
         ylim=rev(range(lc_temp$mag)),
         xlab="Phase",ylab="Magnitude",
         xlim=c(0,2),xaxs='i',cex.axis=1.5,cex.lab=1.5)
    segments(lc_temp$time,
             lc_temp$mag+lc_temp$error,
             lc_temp$time,
             lc_temp$mag-lc_temp$error,col='grey')
    ti <- (1:100)/100
    ti <- c(ti,ti+1)
    m <- PredictAllBand(ti,1,coeffs,tem)
    for(ii in 1:ncol(m)){
        points(ti,m[,ii],type='l',col=bandcol[colnames(m)[ii]])
    }
    legend("bottomleft",names(bandcol),col=bandcol,pch=bandpch,lty=bandpch)
}
