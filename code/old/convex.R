rm(list=ls())
library(ellipse)

n <- 100
X <- cbind(rep(1,n),runif(n,min=1,max=3))
beta <- matrix(c(1,2),ncol=1)
e.sd <- 0.3
y <- X%*%beta + rnorm(n,sd=e.sd)
plot(X[,2],y)
lm.fit <- lm(y~X-1)


A <- t(X)%*%X
C <- -2*t(X)%*%y
beta <- (-1/2)*solve(A,C)
sds <- solve(A)


x <- ellipse(sds,centre=c(beta[1,1],beta[2,1]),level=.99)
ylim <- range(x[,2])
xlim <- range(x[,1])

plot(beta[1,1],beta[2,1],ylim=ylim,xlim=xlim)
points(ellipse(sds,centre=c(beta[1,1],beta[2,1]),level=.95),type='l',col="red",lwd=2)
points(ellipse(sds,centre=c(beta[1,1],beta[2,1]),level=.99),type='l',col="red",lwd=2)
abline(v=0)
abline(h=0)



d <- c(1,0)

sum(d*beta)

z <- solve(A,d)
xstar <- as.numeric(beta - z *(sum(d*beta)/sum(d*z)))
points(xstar[1],xstar[2])
