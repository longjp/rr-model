rm(list=ls())
library(RColorBrewer)
load("template.RData")

templates <- tem$templates
t <- tem$temp_time
ntemp <- nrow(templates)

##plot of templates
ylim <- rev(range(templates))
xlim <- range(t)
pdf("templates.pdf",width=12,height=6)
par(mar=c(5,5,1,1))
plot(0,0,col=0,ylim=ylim,xlim=xlim,xaxs='i',xlab="Phase",ylab="Normalized Magnitude",
     cex.lab=1.3)
for(ii in 1:5){
    points(t,templates[ii,],type='l',lwd=3,col=ii,lty=ii)
}
legend("bottomleft",rownames(templates),col=1:ntemp,lty=1:ntemp,lwd=3,
       cex=1.3)
dev.off()



## ##plot of template derivatives
## dev.new()
## ylim <- range(templatesd)
## xlim <- range(t)
## plot(0,0,col=0,ylim=ylim,xlim=xlim)
## for(ii in 1:5){
##     points(t,templatesd[ii,],type='l',lwd=2)
## }
