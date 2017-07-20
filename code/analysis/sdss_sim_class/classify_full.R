## develop a classifier on well sampled light curves
## which has very good performance
## derive features other than model outputs


rm(list=ls())
set.seed(1234)

## load necessary libraries
library('parallel')
library('multiband')
library('randomForest')
load("../../fit_template/template_sdss.RData")
source("../../fit_template/fit_template.R")
source("../../common/funcs.R")
source("../funcs.R")
source("../../common/plot_funcs.R")

## data source
load("../../data/clean/sdss_sim_class.RData")
load("results.RData")
source("../params.R")



## arguments
##   rho : vector of phases on [0,1) scale
##
## value
## k_hat : 1-step mle of kappa from
##         fisher von mises dist on unit circle
mlefvonmises_simple_k <- function(rho){
    X <- cbind(cos(2*pi*rho),sin(2*pi*rho))
    x_sum <- colSums(X)
    x_sum_norm <- sqrt(sum(x_sum*x_sum))
    Rbar <- x_sum_norm / nrow(X)
    k_hat <- (Rbar*(2-Rbar^2)) / (1-Rbar^2)
    return(k_hat)
}
        




ps <- period_est_FULL[,1] ## just use best fit period


coeffs <- matrix(0,nrow=length(ps),ncol=4)
rss <- rep(0,length(ps))
kappa_feat <- rep(0,length(ps))
phis <- (1:100)/100
for(ii in 1:nrow(coeffs)){
    omega <- 1/ps[ii]
    lc <- TMtoLC(tms_FULL[[ii]])
    rss_phase <- ComputeRSSPhase(lc,omega,tem,phis=phis)
    phi <- phis[which.min(rss_phase)]
    coeffs[ii,] <- ComputeCoeffsPhase(lc,omega,phi,tem)
    pred <- PredictTimeBand(lc[,1],lc[,2],omega,coeffs[ii,],tem)
    rss[ii] <- median(abs(lc[,3] - pred))
    ##rss[ii] <- median((lc[,3] - pred)^2 / (tem$model_error[lc$band]^2 + lc[,4]^2))
    kappa_feat[ii] <- mlefvonmises_simple_k((lc[,1] %% ps[ii])/ps[ii])
    ##rss[ii] <- median(abs(lc[,3] - pred) / sqrt(lc[,4]^2 + 0.03^2))
}
coeffs <- cbind(coeffs,rss,kappa_feat)
colnames(coeffs) <- c("mu","ebv","p2p-gband","phase","rss","kappa")

to_use <- cl=="rr"
summary(rss[to_use])
summary(rss)



### plot kappa versus period
col_pch <- c(1,2)
names(col_pch) <- c("not","rr")
pdf("kappa_period_sdss.pdf")
par(mar=c(5,5,1,1))
plot(ps,kappa_feat,xlab="Period Estimate",ylab="k",cex.lab=1.3,cex.axis=1.3,
     col=col_pch[cl],pch=col_pch[cl])
legend("topleft",c("Not RR","RR"),pch=col_pch,col=col_pch,cex=1.5)
dev.off()





## sort RR Lyrae by RSS, look for any especially bad fits
## looking for cases where the model is failing
rss_mod <- rss
rss_mod[cl!="rr"] <- 0 ## artificially make rss of non-rr 0
ords <- order(rss_mod,decreasing=TRUE)
ii <- 0

ii <- ii + 1
p_est <- ps[ords[ii]]
omega <- 1/p_est
lc <- TMtoLC(tms[[ords[ii]]])
ce2 <- ComputeCoeffs(lc,omega,tem,NN=20)
coeffs[ords[ii],]
plotLC(lc,p_est,coeffs[ords[ii],1:4],tem)
ii
## phase shift issues account for a decent number of actual RR Lyrae with high RSS
## problem is that ComputeCoeffs outputs poor phase estimate, creating large RSS
## possible solution: do grid search on phase after finding period (about 5x longer than current method)
hist(coeffs[ords[1:30],4])



par(mfcol=c(1,3))
plot(ps[to_use],coeffs[to_use,3],ylim=c(0,2))
plot(ps[to_use],coeffs[to_use,5],ylim=c(0,.1))
plot(periods[to_use],ps[to_use],ylim=c(.2,1),xlim=c(.2,1))
abline(a=0,b=1)


## we get 90% of period estimates correct on full light curves, doesn't seem so good
e <- 20 / (60*60*24)
table(abs(ps[to_use]-periods[to_use]) < e)
table(abs(ps[to_use]-periods[to_use]) < e) / sum(to_use)




### find and plot a RR Lyrae with period .5 and a aliased non RR with period .5
tms_FULL_not <- tms_FULL[cl!="rr"]
ps_not <- ps[cl!="rr"]
coeffs_not <- coeffs[cl!="rr",]
ords <- order(abs(ps_not-.5),decreasing=FALSE)
ii <- 2
tm <- tms_FULL_not[[ords[ii]]]
lc <- TMtoLC(tm)
names(lc)[4] <- "error"
p_est <- ps_not[ords[ii]]
pdf("alias_not_rr.pdf",width=12,height=5)
plotLC(lc,p_est,coeffs_not[ords[ii],1:4],tem,
       main=paste0("period estimate=",round(p_est,5),"   kappa=",round(coeffs_not[ords[ii],6],2)))
dev.off()




### find and plot a RR Lyrae with period .5 and a aliased non RR with period .5
tms_FULL_rr <- tms_FULL[cl=="rr"]
ps_rr <- ps[cl=="rr"]
coeffs_rr <- coeffs[cl=="rr",]
ords <- order(abs(ps_rr-.5),decreasing=FALSE)
ii <- 1
tm <- tms_FULL_rr[[ords[ii]]]
lc <- TMtoLC(tm)
names(lc)[4] <- "error"
p_est <- ps_rr[ords[ii]]
pdf("alias_rr_rr.pdf",width=12,height=5)
plotLC(lc,p_est,coeffs_rr[ords[ii],1:4],tem,
       main=paste0("period estimate=",round(p_est,5),"   kappa=",round(coeffs_rr[ords[ii],6],2)))
dev.off()







## get RSS measures, fit RF, get classifier accuracy, example lcs which are misclassified
## possible reasons for missclassification
##   1. rrl with incorrect parameter estimates
##   2. lcs labeled as non-rrl are actually rrl (or look a lot like them)
##   3. need more features


### TODOs

## possible reasons: 1) error distribution of l.c.s with wrong periods has a few wild outliers
##                   2) l.c.s. are poorly sampled
##                   3) code note finding model global min



feats <- cbind(ps,coeffs[,c(3,5,6)])
rf.fit <- randomForest(feats,as.factor(cl))
rf.fit

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


head(preds_not[ords])
ii <- 0



ii <- ii + 1
tm <- tms_FULL_not[[ords[ii]]]
lc <- TMtoLC(tm)
names(lc)[4] <- "error"
p_est <- ps_not[ords[ii]]
plotLC(lc,p_est,coeffs_not[ords[ii],1:4],tem,
       main=paste0("prob RR=",round(preds_not[ords[ii]],3),"  period=",round(p_est,3)))
ii

### summary of non--rrlyrae classified as RR Lyrae
## 3 likely rr lyrae misclassified as not RR
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


ii <- 0



ii <- ii + 1
tm <- tms_FULL_rr[[ords[ii]]]
lc <- TMtoLC(tm)
names(lc)[4] <- "error"
p_est <- ps_rr[ords[ii]]
plotLC(lc,p_est,coeffs_rr[ords[ii],1:4],tem,
       main=paste0("prob RR=",round(preds_rr[ords[ii]],3),"  period=",round(p_est,3)))


