## develop a classifier on well sampled light curves
## which has very good performance
## derive features other than model outputs


rm(list=ls())
set.seed(1234)

## load necessary libraries
library('parallel')
library('multiband')
library('randomForest')
load("../../fit_template/template.RData")
source("../../fit_template/fit_template.R")
source("../../common/funcs.R")
source("../funcs.R")


## data source
load("../../data/clean/sdss_sim_class.RData")
load("results.RData")
source("../params.R")








ps <- period_est_FULL[,1] ## just use best fit period


coeffs <- matrix(0,nrow=length(ps),ncol=4)
rss <- rep(0,length(ps))
for(ii in 1:nrow(coeffs)){
    omega <- 1/ps[ii]
    lc <- TMtoLC(tms_FULL[[ii]])
    coeffs[ii,] <- ComputeCoeffs(lc,omega,tem,NN=20,use.errors=FALSE)
    pred <- PredictTimeBand(lc[,1],lc[,2],omega,coeffs[ii,],tem)
    rss[ii] <- median(abs(lc[,3] - pred))
}
coeffs <- cbind(coeffs,rss)
colnames(coeffs) <- c("mu","ebv","p2p-gband","phase","rss")

to_use <- cl=="rr"
summary(rss[to_use])



par(mfcol=c(1,3))
plot(ps[to_use],coeffs[to_use,3],ylim=c(0,2))
plot(ps[to_use],coeffs[to_use,5],ylim=c(0,.1))
plot(periods[to_use],ps[to_use],ylim=c(.2,1),xlim=c(.2,1))
abline(a=0,b=1)


## we get 90% of period estimates correct on full light curves, doesn't seem so good
e <- 20 / (60*60*24)
table(abs(ps[to_use]-periods[to_use]) < e)


## get RSS measures, fit RF, get classifier accuracy, example lcs which are misclassified
## possible reasons for missclassification
##   1. rrl with incorrect parameter estimates
##   2. lcs labeled as non-rrl are actually rrl (or look a lot like them)
##   3. need more features


### TODOs

## possible reasons: 1) error distribution of l.c.s with wrong periods has a few wild outliers
##                   2) l.c.s. are poorly sampled
##                   3) code note finding model global min



feats <- cbind(ps,coeffs[,c(3,5)])
rf.fit <- randomForest(feats,as.factor(cl))


preds <- predict(rf.fit,type='prob')[,2]
preds_cl <- predict(rf.fit)



ds <- lapply(unique(cl),function(x){density(preds[cl==x])})
names(ds) <- unique(cl)
xlim <- range(unlist(lapply(ds,function(d){d$x})))
ylim <- range(unlist(lapply(ds,function(d){d$y})))
plot(0,0,xlim=xlim,ylim=ylim,col=0)
for(ii in 1:length(ds)){
    points(ds[[ii]]$x,ds[[ii]]$y,type='l',col=ii)
}



## order non--rrlyrae from most likely to be rr to least likely
## then plot
tms_FULL_not <- tms_FULL[cl!="rr"]
preds_not <- preds[cl!="rr"]
ps_not <- ps[cl!="rr"]
coeffs_not <- coeffs[cl!="rr",]
ords <- order(preds_not,decreasing=TRUE)


preds_not[ords]
ii <- 0



ii <- ii + 1
tm <- tms_FULL_not[[ords[ii]]]
lc <- TMtoLC(tm)
names(lc)[4] <- "error"
p_est <- ps_not[ords[ii]]




colpch <- 1:5
names(colpch) <- names(tem$betas)

lc1 <- lc
lc1[,1] <- (lc$time %% p_est)/p_est
lc2 <- lc1
lc2[,1] <- lc1[,1] + 1
lc_temp <-rbind(lc1,lc2)
plot(lc_temp$time,lc_temp$mag,
     col=colpch[lc_temp$band],pch=colpch[lc_temp$band],
     ylim=rev(range(lc_temp$mag)),
     xlab="time",ylab="magnitude",
     xlim=c(0,2),xaxs='i',main=paste0("prob RR=",round(preds_not[ords[ii]],3),"  period=",round(p_est,3)))
segments(lc_temp$time,
         lc_temp$mag+lc_temp$error,
         lc_temp$time,
         lc_temp$mag-lc_temp$error)
ti <- (1:100)/100
ti <- c(ti,ti+1)
m <- PredictAllBand(ti,1,coeffs_not[ords[ii],1:4],tem)
for(jj in 1:length(tem$betas)){
    points(ti,m[,jj],type='l',col=colpch[names(tem$betas)[jj]])
}
head(m)
summary(m)


### summary of non--rrlyrae classified as RR Lyrae
## 2 likely rr lyrae misclassified as not RR
## several (maybe 10) with aliasing at .5 or .333 days
## a few periodic variables that are not RRab (maybe c, or other classes)
## few ambiguous cases

##### conclusion: find method for addressing aliasing, mostly occurs at .5 and .333, but some at 1, .25

##### aliasing occurs because of observing star at the same time on every night
### - if there is some other pattern of observing, will have different aliasing
### - check if DES is same time each night, how can we address this







### examine RRL classified as not
## order non--rrlyrae from most likely to be rr to least likely
## then plot
tms_FULL_rr <- tms_FULL[cl=="rr"]
preds_rr <- preds[cl=="rr"]
ps_rr <- ps[cl=="rr"]
coeffs_rr <- coeffs[cl=="rr",]
ords <- order(preds_rr)


preds_rr[ords]
ii <- 0



ii <- ii + 1
tm <- tms_FULL_rr[[ords[ii]]]
lc <- TMtoLC(tm)
names(lc)[4] <- "error"
p_est <- ps_rr[ords[ii]]

colpch <- 1:5
names(colpch) <- names(tem$betas)

lc1 <- lc
lc1[,1] <- (lc$time %% p_est)/p_est
lc2 <- lc1
lc2[,1] <- lc1[,1] + 1
lc_temp <-rbind(lc1,lc2)
plot(lc_temp$time,lc_temp$mag,
     col=colpch[lc_temp$band],pch=colpch[lc_temp$band],
     ylim=rev(range(lc_temp$mag)),
     xlab="time",ylab="magnitude",
     xlim=c(0,2),xaxs='i',main=paste0("prob RR=",round(preds_rr[ords[ii]],3),"  period=",round(p_est,3)))
segments(lc_temp$time,
         lc_temp$mag+lc_temp$error,
         lc_temp$time,
         lc_temp$mag-lc_temp$error)
ti <- (1:100)/100
ti <- c(ti,ti+1)
m <- PredictAllBand(ti,1,coeffs_rr[ords[ii],1:4],tem)
for(jj in 1:length(tem$betas)){
    points(ti,m[,jj],type='l',col=colpch[names(tem$betas)[jj]])
}
head(m)
summary(m)



## 2 with long periods (but looked pretty good)
## 4 insane outliers


## no rrab have periods < .4, so grid search over .4 to 1 would be faster and reduce aliasing






plot(density(periods[cl=="rr"],bw="SJ"))

plot(density(periods[cl=="rr"],bw=.001))
a

plot(density(period_est_FULL))


sum(period_est_FULL[,1] < .52 & period_est_FULL[,1] > .48)


to_use <- (period_est_FULL[,1] < .501) & (period_est_FULL[,1] > .499)
table(cl[to_use],preds_cl[to_use])


sum(period_est_FULL[,1] < .51 & period_est_FULL[,1] > .49)
sum(period_est_FULL[,1] < .51 & period_est_FULL[,1] > .49 & preds_cl!="rr")
sum(period_est_FULL[,1] < .51 & period_est_FULL[,1] > .49 & cl!="rr")


## is star RR or not
## is star sampled 1 day exact cadence
## if RR, does star have period .5


## if star is NOT an RR and regularly sampled
## - when folded at .5 days, all same phase space, (amp, phase) become "random", very large uncertainties
## rss = photometric error (this is assuming mu and ebv can model colors)
## - when folded NOT at .5 days, photometric measurements well spread in phase space, amplitude near 0, phase random
##   rss = photometric error (this is assuming mu and ebv can model colors)




## if star IS and RR and regularly sampled
## - when folded at .5 days, all same phase space, (amp, phase) become "random", very large uncertainties
## rss = photometric error + scatter due to variation (this is assuming mu and ebv can model colors), unless star actually
##       has period near .5 days
## - when not .5 days, or other alias, will look poor unless period is correct

## AN RR LYRAE WITH PERIOD .5 DAYS AND REGULARLY SAMPLED IS IDENTICAL TO CONSTANT STAR WITH RR CONSISTENT COLORS
## there are many RRL with periods near .5 days that look very nice when folded



### PLOT RR LYRAE with period estimates (or true periods) near 0.5





### examine RRL classified as not
## order non--rrlyrae from most likely to be rr to least likely
## then plot
tms_FULL_rr <- tms_FULL[cl=="rr"]
preds_rr <- preds[cl=="rr"]
ps_rr <- ps[cl=="rr"]
coeffs_rr <- coeffs[cl=="rr",]
ords <- order(abs(ps_rr - .5))


preds_rr[ords]
ii <- 0



ii <- ii + 1
tm <- tms_FULL_rr[[ords[ii]]]
lc <- TMtoLC(tm)
names(lc)[4] <- "error"
p_est <- ps_rr[ords[ii]]

colpch <- 1:5
names(colpch) <- names(tem$betas)

lc1 <- lc
lc1[,1] <- (lc$time %% p_est)/p_est
lc2 <- lc1
lc2[,1] <- lc1[,1] + 1
lc_temp <-rbind(lc1,lc2)
plot(lc_temp$time,lc_temp$mag,
     col=colpch[lc_temp$band],pch=colpch[lc_temp$band],
     ylim=rev(range(lc_temp$mag)),
     xlab="time",ylab="magnitude",
     xlim=c(0,2),xaxs='i',main=paste0("prob RR=",round(preds_rr[ords[ii]],3),"  period=",round(p_est,3)))
segments(lc_temp$time,
         lc_temp$mag+lc_temp$error,
         lc_temp$time,
         lc_temp$mag-lc_temp$error)
ti <- (1:100)/100
ti <- c(ti,ti+1)
m <- PredictAllBand(ti,1,coeffs_rr[ords[ii],1:4],tem)
for(jj in 1:length(tem$betas)){
    points(ti,m[,jj],type='l',col=colpch[names(tem$betas)[jj]])
}
head(m)
summary(m)





