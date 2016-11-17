rm(list=ls())
library(RColorBrewer)
load("template.RData")

templates <- tem$templates
t <- tem$temp_time

##plot of templates
ylim <- range(templates)
xlim <- range(t)
pdf("templates2.pdf",width=8,height=4)
plot(0,0,col=0,ylim=ylim,xlim=xlim)
for(ii in 1:5){
    points(t,templates[ii,],type='l',lwd=2)
}
dev.off()

##plot of template derivatives
dev.new()
ylim <- range(templatesd)
xlim <- range(t)
plot(0,0,col=0,ylim=ylim,xlim=xlim)
for(ii in 1:5){
    points(t,templatesd[ii,],type='l',lwd=2)
}
