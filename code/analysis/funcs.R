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
    ps <- SeparateBest(ps,0.0002,topN)
    return(ps)
}

FitLombParallel <- function(ii,tms,omegas,topN){
    print(ii)
    tm <- tms[[ii]]
    p_grid <- rev(1/omegas)
    out <- pgls(tm,periods=p_grid,BCD_flag=FALSE)
    rss <- out$rss_ls
    ps <- p_grid[sort_local_min(1:length(rss),rss)]
    ps <- SeparateBest(ps,0.0002,topN)
    return(ps)
}
