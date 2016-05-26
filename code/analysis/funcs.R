FitTemplateParallel <- function(ii,tms,omegas,tem,NN){
    print(ii)
    rss <- FitTemplate(tms[[ii]],omegas,tem,NN=NN)
    ps <- 1/omegas[sort_local_min(1:length(rss),rss)]
    ps <- SeparateBest(ps,0.0002,5)
    return(ps)
}


GetFreqs <-function(period_min,period_max,freq_del = 0.1/4000){
    freq_max <- 1/period_min
    freq_min <- 1/period_max
    return(seq(freq_min, freq_max, freq_del))
}
