rm(list=ls())
load("features.RData")
features$class <- as.factor(features$class)
fnames <- rownames(features)

lum <- rep(0,nrow(features))
for(ii in 1:nrow(features)){
    lum[ii] <- mean(read.table(fnames[ii])[,2])
}


plot(log(features$p,base=10),lum,ylim=c(max(lum),min(lum)),
     col=c("orange", "blue", "black","green3","red")[features$class],
     pch=(1:5)[features[,5]],xlab="log(period)",
     ylab="Luminosity (mean I band magnitude)")



n <- 100

x1 <- rexp(n,rate=1/50) + rnorm(n,mean=130,sd=10)
x2 <- rexp(n,rate=1/50) + rnorm(n,mean=140,sd=10)
x3 <- rexp(n,rate=1/50) + rnorm(n,mean=150,sd=10)

a <- 2
lum1 <- -a*log(x1,base=10) + 20
lum2 <- -a*log(x2,base=10) + 19
lum3 <- -a*log(x3,base=10) + 18

dat <- as.data.frame(cbind(period=c(x1,x2,x3),luminosity=c(lum1,lum2,lum3)))
plot(log(dat$period,base=10),dat$luminosity)
