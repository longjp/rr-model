LCtoTM <- function(lc){
    lc[,1] <- lc[,1] - min(lc[,1])
    levs <- as.character(levels(lc$b))
    levs <- levs[order(levs)]
    tm <- list()
    for(ii in 1:length(levs)){
        tm[[ii]] <- lc[lc$b == levs[ii],c("time","mag","sigma")]
        names(tm[[ii]]) <- c("time","mag","sigma")
    }
    names(tm) <- levs
    return(tm)
}

TMtoLC <- function(tm){
    nb <- vapply(tm,nrow,c(0))
    bs <- rep.int(names(tm),nb)
    lc <- do.call(rbind,tm)
    lc <- data.frame(time=lc[,1],band=bs,mag=lc[,2],sigma=lc[,3])
    return(lc)
}

    
