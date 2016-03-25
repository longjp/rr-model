SplitData <- function(tm){
    levs <- levels(tm$band)
    lc <- list()
    for(ii in 1:length(levs)){
        lc[[ii]] <- tm[tm$band==levs[ii],c("time","mag","error")]
    }
    names(lc) <- levs
    return(lc)
}
