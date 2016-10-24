rm(list=ls())
set.seed(1234)
load("features.RData")
features$p <- log(features$p,base=10)
names(features)[1:3] <- c("log(period)","skew","amp")
features[,5] <- NULL


pdf("figs/pairs_cluster.pdf",width=6,height=4)
par(mar=c(rep(2,4)))
pairs(features[1:4])
dev.off()



features <- as.matrix(features)
features.c <- scale(features)
f.d <- dist(features.c)
fit <- hclust(f.d)



pdf("figs/dend.pdf",width=8,height=4)
par(mar=c(1,1,0,1))
plot(fit,labels=FALSE,hang=0,main="",axes=FALSE)
dev.off()


for(ii in 2:5){
    pdf(paste("figs/dend",ii,".pdf",sep=""),width=8,height=4)
    par(mar=c(1,1,0,1))
    plot(fit,labels=FALSE,hang=0,main="",axes=FALSE)
    rect.hclust(fit, k=ii, border="red")
    dev.off()
    pdf(paste("figs/pairs_cluster",ii,".pdf",sep=""),width=6,height=4)
    fit.cut <- cutree(fit,k=ii)
    par(mar=c(rep(2,4)))
    pairs(features,col=fit.cut,pch=fit.cut)
    dev.off()
}




