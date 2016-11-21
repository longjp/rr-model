rm(list=ls())
library(RColorBrewer)
load("template.RData")

templates <- tem$templates
t <- tem$temp_time
ntemp <- nrow(templates)


band_mark <- 1:6
names(band_mark) <- c("u","g","r","i","z","Y")


##plot of templates
ylim <- rev(range(templates))
xlim <- range(t)
pdf("templates.pdf",width=12,height=6)
par(mar=c(5,5,1,1))
plot(0,0,col=0,ylim=ylim,xlim=xlim,xaxs='i',xlab="Phase",ylab="Normalized Magnitude",
     cex.lab=2)
for(ii in 1:nrow(templates)){
    points(t,templates[ii,],type='l',lwd=3,
           col=band_mark[rownames(templates)[ii]],
           lty=band_mark[rownames(templates)[ii]])
}
to_use <- names(band_mark) %in% rownames(templates)
legend("bottomleft",names(band_mark)[to_use],
       col=band_mark[to_use],lty=band_mark[to_use],lwd=3,cex=2)
dev.off()



## ##plot of template derivatives
## dev.new()
## ylim <- range(templatesd)
## xlim <- range(t)
## plot(0,0,col=0,ylim=ylim,xlim=xlim)
## for(ii in 1:5){
##     points(t,templatesd[ii,],type='l',lwd=2)
## }
