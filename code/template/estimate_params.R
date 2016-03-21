rm(list=ls())
load("sim.RData")
load("make_template.RData")

temp_time <- seq(0,1,length.out=ncol(templates))




AugmentData <- function(tm,dust,betas,use.errors=FALSE){
    tm <- tm[order(tm[,2]),]
    nb <- table(tm[,2])
    tm$dust <- rep.int(dust,nb)
    tm$mag <- tm$mag - rep.int(betas,nb)
    tm$band <- NULL
    if(!use.errors){
        tm$error <- 1
    }
    return(list(tm=tm,nb=nb))
}


ii <- 1
dat <- AugmentData(tms[[ii]],dust,betas)


par(mfcol=c(2,1))
plot((dat[[1]]$time %% period[ii]) / period[ii],dat[[1]]$mag - dat[[1]]$dust*d[ii])
plot((tms[[ii]]$time %% period[ii]) / period[ii],tms[[1]]$mag)





ConstructGamma <- function(t,nb,phi,per,templates,temp_time){
    t <- (t/per + phi) %% 1
    bix <- c(0,cumsum(nb))
    gammaf <- rep(0,length(t))
    for(jj in 1:(length(bix)-1)){
        ix1 <- (bix[jj]+1)
        ix2 <- bix[jj+1]
        gammaf[ix1:ix2] <-  approx(temp_time,templates[jj,],xout=t[ix1:ix2])$y
    }
    return(gammaf)
}


## computes beta
ComputeBeta <- function(m,dust,gammaf){
    X <- cbind(1,dust,gammaf)
    B <- t(X)%*%X
    d <- t(X)%*%m
    z <- solve(B,d)
   return(z)
}



t <- dat[[1]]$time
nb <- dat[[2]]
phi <- phase[ii]
per <- period[ii]
m <- dat[[1]]$mag

gammaf <- ConstructGamma(t,nb,phi,per,templates,temp_time)


pred <- gammaf*a[ii] + d[ii]*dat[[1]]$dust + rep(alpha[ii],length(gammaf)) + rep.int(betas,nb)

tms[[ii]] <- tms[[ii]][order(tms[[ii]][,2]),]

plot((tms[[ii]][,1] %% per) / per,tms[[ii]][,3])
points((tms[[ii]][,1] %% per) / per,pred,col='red')




params <- ComputeBeta(m,dat[[1]]$dust,gammaf)
