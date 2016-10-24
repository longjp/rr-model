rm(list=ls())
set.seed(1234)
library(tree)
library(tables)
library(MASS)
source("plot-partition.R")

generate_data <- function(n){
    x <- runif(n)
    y <- runif(n)
    label <- rep(1,n)
    label <- label + 1*(y < 2*x)
    label <- label + 1*((y<0.4) & (y < 2*x))
    replace.original <- rbinom(n,1,.1)
    label <- (1 - replace.original)*label + replace.original*sample(1:3,n,replace=TRUE)
    dat <- data.frame(class=c("black","orange","blue")[label],x=x,y=y)
    return(dat)
}


## generate_data <- function(n){
##     x <- runif(n)
##     y <- runif(n)
##     label <- rep(1,n)
##     label <- label + 1*(y < 2*x)
##     label <- label + 1*((y<0.4) & (y < 2*x))
##     replace.original <- rbinom(n,1,.1)
##     label <- (1 - replace.original)*label + replace.original*sample(1:3,n,replace=TRUE)
##     dat <- data.frame(class=c("black","orange","blue")[label],x=x,y=y)
##     return(dat)
## }

n <- 1000
dat <- generate_data(n)
test <- generate_data(100)


############# EXPERIMENT BEGINS HERE

head(dat)
lims <- c(range(dat$x),range(dat$y))

h <- .2
plot(0,0,xlim=c(-.5,1.5),ylim=c(-.5,1.5),col=0,xlab="Feature x",ylab="Feature y",cex.lab=1.5)
points(dat$x,dat$y,col=as.character(dat$class),pch=as.numeric(dat$class))
for(jj in unique(dat$class)){
    to_use <- dat$class==jj
    fit.kde <- kde2d(dat$x[to_use],dat$y[to_use],n=50,h=h,lims=c(-.5,1.5,-.5,1.5))
    contour(fit.kde,levels  =  c(0.05, 0.1, 0.2, 0.4),add=TRUE,col=jj)
}



## for(jj in unique(dat$class)){
##     to_use <- dat$class==jj
##     fit.kde <- kde2d(dat$x[to_use],dat$y[to_use],n=50,h=h,lims=c(-.5,1.5,-.5,1.5))
##     contour(fit.kde,levels  =  c(0.05, 0.1, 0.2, 0.4),add=TRUE,col=jj)
## }




############## EXPERIMENT ENDS HERE















fit.tree <- tree(class~.,data=dat)






pdf("figs/sim_features.pdf",width=6,height=6)
par(mar=c(4.5,4.5,1,1))
plot(dat$x,dat$y,col=as.character(dat$class),xlab="Feature x",ylab="Feature y",cex.lab=1.5,
     pch=as.numeric(dat$class))
dev.off()


## plot sequence of trees
n.best <- max(prune.tree(fit.tree)$size)
for(ii in 2:n.best){
    fit <- prune.tree(fit.tree,best=ii)
    pdf(paste("figs/tree_",ii,".pdf",sep=""),width=12,height=6)
    par(mfcol=c(1,2),mar=c(4.5,4.5,1,1))
    plot(dat$x,dat$y,col=as.character(dat$class),xlab="Feature x",ylab="Feature y",cex.lab=1.5,
         pch=as.numeric(dat$class))
    partition.tree(fit,cex=1.5,add=TRUE,ordvars=c("x","y"),font=2)
    plot(fit,cex.lab=2,type="uniform")
    text(fit,digits=1,cex=1.5,font=2)
    dev.off()
}



tblr <- tabular((Truth=test$class) ~ (Predicted = predict(fit.tree,newdata=test,type="class")))
sink("figs/confusion_sim.tex")
latex(tblr)
sink()
