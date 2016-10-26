rm(list=ls())


## generate bivariate normal
p <- 2/3

rnormal_mix <- function(n=1000){
    m <- c(4,0)
    s <- c(1,1)
    ind <- rbinom(n,1,p) + 1
    return(m[ind] + s[ind]*rnorm(n))
}

normal_mix <- function(x){
    return(p*dnorm(x,mean=0,sd=1) + (1-p)*dnorm(x,mean=4,sd=1))
}



ylim <- c(0,max(normal_mix(seq(from=-4,to=8,length.out=1000))))
xlim <- c(-4,8)
x <- rnormal_mix()
sj <- bw.SJ(x)
bws <- c(sj/10,sj,5*sj)
sj <- bw.SJ(x)

for(jj in 1:2){
    x <- rnormal_mix()
    for(ii in 1:length(bws)){
        pdf(paste0("density_figs/density",ii,"_",jj,".pdf"))
        par(mar=c(5,5,1,1))
        plot(density(x,bw=bws[ii]),main="",col='red',lty=2,
             xlab=paste0("x (bandwidth=",round(bws[ii],2),")"),
             ylim=ylim,xlim=xlim,lwd=2)
        curve(normal_mix,add=TRUE,lwd=2)
        points(x,rep(0,length(x)),col="#00000030")
        dev.off()
    }
}

### smallest bandwidth, but with larger sample size
x <- rnormal_mix(n=10000)
pdf(paste0("density_figs/density",1,"_largen.pdf"))
par(mar=c(5,5,1,1))
plot(density(x,bw=bws[1]),main="",col='red',lty=2,
     xlab=paste0("x (bandwidth=",round(bws[1],2),")"),
     ylim=ylim,xlim=xlim,lwd=2)
curve(normal_mix,add=TRUE,lwd=2)
points(x,rep(0,length(x)),col="#00000030")
dev.off()


file.copy(from=list.files('density_figs',full.names=TRUE),
          to=paste0('../figs/',list.files('density_figs')),
          overwrite=TRUE)

