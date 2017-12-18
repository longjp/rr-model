## returns dust corrected lightcurve
##    arguments:
##          tm : lightcurve
##         ebv : dust for lightcurve tm
##         tem : templates, contains extinction law in tem$dust
DustCorrect <- function(tm,ebv,tem){
    if(mean(names(tm) %in% names(tem$dust))!=1){
        print("lightcurve has bands: ",names(tm))
        print("template only has extinction law for:",names(tem$dust))
        stop()
    }
    for(jj in 1:length(tm)){
        tm[[jj]][,2] <- tm[[jj]][,2] - tem$dust[names(tm)[jj]]*ebv
    }
    return(tm)
}

GetFreqs <-function(period_min,period_max,freq_del = 0.1/4000){
    freq_max <- 1/period_min
    freq_min <- 1/period_max
    return(seq(freq_min, freq_max, freq_del))
}

FitTemplateParallel <- function(ii,tms,omegas,tem,NN=5,use.errors=TRUE,use.dust=TRUE,topN=topN){
    print(ii)
    lc <- TMtoLC(tms[[ii]])
    rss <- FitTemplate(lc,omegas,tem,NN=NN,use.errors=use.errors,use.dust=use.dust)
    ps <- 1/omegas[sort_local_min(1:length(rss),rss)]
    ps <- SeparateBest(ps,0.01,topN)
    return(ps)
}

FitLombParallel <- function(ii,tms,omegas,topN){
    print(ii)
    tm <- tms[[ii]]
    ## don't use errors with lomb because model is crude approximation
    for(jj in 1:length(tm)){
        tm[[jj]][,3] <- 1
    }
    p_grid <- rev(1/omegas)
    out <- pgls(tm,periods=p_grid,BCD_flag=FALSE)
    rss <- out$rss_ls
    ps <- p_grid[sort_local_min(1:length(rss),rss)]
    ps <- SeparateBest(ps,0.0002,topN)
    return(ps)
}
