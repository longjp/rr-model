partition.tree <- function (tree, label = "yval", add = FALSE, ordvars, ...){
    ptXlines <- function(x, v, xrange, xcoord = NULL, ycoord = NULL, 
        tvar, i = 1L) {
        if (v[i] == "<leaf>") {
            y1 <- (xrange[1L] + xrange[3L])/2
            y2 <- (xrange[2L] + xrange[4L])/2
            return(list(xcoord = xcoord, ycoord = c(ycoord, y1, 
                y2), i = i))
        }
        if (v[i] == tvar[1L]) {
            xcoord <- c(xcoord, x[i], xrange[2L], x[i], xrange[4L])
            xr <- xrange
            xr[3L] <- x[i]
            ll2 <- Recall(x, v, xr, xcoord, ycoord, tvar, i + 
                1L)
            xr <- xrange
            xr[1L] <- x[i]
            return(Recall(x, v, xr, ll2$xcoord, ll2$ycoord, tvar, 
                ll2$i + 1L))
        }
        else if (v[i] == tvar[2L]) {
            xcoord <- c(xcoord, xrange[1L], x[i], xrange[3L], 
                x[i])
            xr <- xrange
            xr[4L] <- x[i]
            ll2 <- Recall(x, v, xr, xcoord, ycoord, tvar, i + 
                1L)
            xr <- xrange
            xr[2L] <- x[i]
            return(Recall(x, v, xr, ll2$xcoord, ll2$ycoord, tvar, 
                ll2$i + 1L))
        }
        else stop("wrong variable numbers in tree.")
    }
    if (inherits(tree, "singlenode")) 
        stop("cannot plot singlenode tree")
    if (!inherits(tree, "tree")) 
        stop("not legitimate tree")
    frame <- tree$frame
    leaves <- frame$var == "<leaf>"
    var <- ordvars
    if (length(var) > 2L || length(var) < 1L) 
        stop("tree can only have one or two predictors")
    nlevels <- sapply(attr(tree, "xlevels"), length)
    if (any(nlevels[var] > 0L)) 
        stop("tree can only have continuous predictors")
    x <- rep(NA, length(leaves))
    x[!leaves] <- as.double(substring(frame$splits[!leaves, "cutleft"], 
        2L, 100L))
    m <- model.frame(tree)
    if (!missing(ordvars)) {
        ind <- match(var, ordvars)
        if (any(is.na(ind))) 
            stop("unmatched names in vars")
        var <- ordvars[sort(ind)]
    }
    lab <- frame$yval[leaves]
    if (is.null(frame$yprob)) 
        lab <- format(signif(lab, 3L))
    else if (match(label, attr(tree, "ylevels"), nomatch = 0L)) 
        lab <- format(signif(frame$yprob[leaves, label], 
                             3L))
    rx <- range(m[[var[1L]]])
    rx <- rx + c(-0.025, 0.025) * diff(rx)
    rz <- range(m[[var[2L]]])
    ##rz <- c(0,1)
    rz <- rz + c(-0.025, 0.025) * diff(rz)
    xrange <- c(rx, rz)[c(1, 3, 2, 4)]
    xcoord <- NULL
    ycoord <- NULL
    xy <- ptXlines(x, frame$var, xrange, xcoord, ycoord, 
                   var)
    xx <- matrix(xy$xcoord, nrow = 4L)
    yy <- matrix(xy$ycoord, nrow = 2L)
    if (!add) 
        plot(rx, rz, xlab = var[1L], ylab = var[2L], type = "n", 
             xaxs = "i", yaxs = "i", ...)
    segments(xx[1L, ], xx[2L, ], xx[3L, ], xx[4L, ])
    text(yy[1L, ], yy[2L, ], as.character(lab), ...)
}
