rm(list=ls())
set.seed(1234)
library(rpart)
library(rpart.plot)
library(xtable)
library(tables)
load('features.RData')


features$class <- as.factor(features$class)
rownames(features) <- 1:nrow(features)
head(features)

names(features)
features$p <- log(features$p,base=10)
names(features)[1:3] <- c("log(period)","skew","amp")

head(features)


pdf("figs/pairs.pdf",width=6,height=4)
pairs(features[1:4],
      col = c("orange", "blue", "black","green3","red")[features[,5]],
      oma=c(2,2,2,12),pch=(1:5)[features[,5]])
# allow plotting of the legend outside the figure region 
# (ie within the space left by making the margins big)
par(xpd=TRUE)
legend(0.83, 0.75, as.vector(unique(features[,5])),  
       col=c("orange", "blue", "black","green3","red"),
       pch=(1:5),cex=1)
dev.off()



## convert period back to non-log
features[,1] <- 10^features[,1]
names(features)[1] <- "period"

## randomize order
features <- features[sample(1:nrow(features),nrow(features)),]
n.train <- 400
features.train <- features[1:n.train,]
features.test <- features[-(1:n.train),]

head(features.train)


fit.rpart <- rpart(class ~ .,data=features.train)

pdf("figs/tree.pdf",width=10,height=5)
par(mar=c(0,1,0,1))
prp(fit.rpart,extra=2,compress=FALSE,varlen=0)
dev.off()


## print confusion matrix
tblr <- tabular((Truth=features.test$class) ~ (Predicted = predict(fit.rpart,features.test,type="class")))
sink("figs/confusion.tex")
latex(tblr)
sink()



plot(log(features$amp,base=10),log(features$p2p_scatter,base=10),col=as.numeric(as.factor(features[,5])))



