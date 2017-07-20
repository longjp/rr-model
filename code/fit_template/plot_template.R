rm(list=ls())
library(RColorBrewer)
load("template_sdss.RData")

templates <- tem$templates
t <- tem$temp_time
ntemp <- nrow(templates)


bandmark <- 1:6
names(bandmark) <- c("u","g","r","i","z","Y")
bandcol <- c("dodgerblue3","green","red",
             "mediumorchid1","black","peachpuff4")
names(bandcol) <- c("u","g","r","i","z","Y")

##plot of templates
ylim <- rev(range(templates))
xlim <- range(t)
pdf("templates.pdf",width=12,height=6)
par(mar=c(5,5,1,1))
plot(0,0,col=0,ylim=ylim,xlim=xlim,xaxs='i',xlab="Phase",ylab="Normalized Magnitude",
     cex.lab=2)
for(ii in 1:nrow(templates)){
    points(t,templates[ii,],type='l',lwd=3,
           col=bandcol[rownames(templates)[ii]],
           lty=bandmark[rownames(templates)[ii]])
}
to_use <- names(bandmark) %in% rownames(templates)
legend("bottomleft",names(bandmark)[to_use],
       col=bandcol[to_use],lty=bandmark[to_use],lwd=3,cex=2)
dev.off()



## ##plot of template derivatives
## dev.new()
## ylim <- range(templatesd)
## xlim <- range(t)
## plot(0,0,col=0,ylim=ylim,xlim=xlim)
## for(ii in 1:5){
##     points(t,templatesd[ii,],type='l',lwd=2)
## }
