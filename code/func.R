get_freqs <-function(period_min,period_max,freq_del = 0.1/4000){
    freq_max <- 1/period_min
    freq_min <- 1/period_max
    return(2 * pi * seq(freq_min, freq_max, freq_del))
}

is_local_min <- function(x){
    x.max <- max(x)
    x.extend <- c(x.max+1,x,x.max+1)
    x.right <- x - x.extend[-c(1,2)]
    x.left <- x - x.extend[1:length(x)]
    return(x.right < 0 & x.left < 0)
}

sort_local_min <- function(x,y){
    loc <- is_local_min(y)
    x <- x[loc]
    y <- y[loc]
    ords <- order(y)
    x <- x[ords]
    return(x)
}
