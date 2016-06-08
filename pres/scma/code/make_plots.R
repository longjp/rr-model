rm(list=ls())
library(multiband)
source("../../../code/common/funcs.R")

cat <- read.table("stripe82candidateVar_v1.1.dat",header=TRUE)

make_plot <- function(lc,period,id=NULL){
    names(lc) <- c("time","band","mag","sigma")
    bands <- unique(lc$band)
    bands <- as.character(bands)
    bands <- sort(bands)

    ylim <- range(lc$mag)
    xlim <- c(0,1)
    pdf(paste0("../figs/folded_",id,".pdf"),width=12,height=6)
    par(mar=c(5,5,1,1))
    plot(0,0,col=0,ylim=rev(ylim),xlim=xlim,xlab=paste0("Phase (period = ",round(period,2),")"),
         ylab="Magnitude",cex.lab=1.5,xaxs='i')
    for(jj in 1:length(bands)){
        temp <- lc[lc$band==bands[jj],]
        segments((temp$time %% period)/period,temp$mag-temp$sigma,
        (temp$time %% period)/period,temp$mag+temp$sigma,col='grey')
        points((temp$time %% period)/period,temp$mag,col=jj,pch=jj,cex=1.5)
    }
    legend("bottomleft",paste0(bands," Band"),col=1:length(bands),pch=1:length(bands),cex=1.5)
    dev.off()
    
    ylim <- range(lc$mag)
    xlim <- range(lc$time)
    pdf(paste0("../figs/unfolded_",id,".pdf"),width=12,height=6)
    par(mar=c(5,5,1,1))
    plot(0,0,col=0,ylim=rev(ylim),xlim=xlim,xlab="Time",ylab="Magnitude",cex.lab=1.5,xaxs='i')
    for(jj in 1:length(bands)){
        temp <- lc[lc$band==bands[jj],]
        segments(temp$time,temp$mag-temp$sigma,temp$time,temp$mag+temp$sigma,col='grey')
        points(temp$time,temp$mag,col=jj,pch=jj,cex=1.5)
    }
    legend("bottomleft",paste0(bands," Band"),col=1:length(bands),pch=1:length(bands),cex=1.5)
    dev.off()
}


rrlyrae <- read.table("apj326724t2_mrt.txt",skip=42)


id <- 4183016
lc <- read.table(paste0("AllLCs/LC_",id,".dat"))
period <- cat$P[cat$ID==id]
make_plot(lc,period,id)

id <- 4099
lc <- read.table(paste0("AllLCs/LC_",id,".dat"))
period <- rrlyrae[rrlyrae$V1==id,3]
make_plot(lc,period,id)


id <- 7904669
lc <- read.table(paste0("AllLCs/LC_",id,".dat"))
period <- cat$P[cat$ID==id]
make_plot(lc,period,id)

id <- 13350
lc <- read.table(paste0("AllLCs/LC_",id,".dat"))
period <- rrlyrae[rrlyrae$V1==id,3]
make_plot(lc,period,id)

lc <- lc[sample(1:304,20),]
make_plot(lc,period,paste0(id,"down"))



lc <- read.table("PS1_sample_LCs/ps1_-3266288406308450226cleaned.txt",
                 header=TRUE)
lc <- lc[,c(3,2,4,1)]
names(lc) <- c("time","band","mag","sigma")
tm <- LCtoTM(lc)
out <- pgls(tm,period_min=0.2,period_max=1,BCD_flag=FALSE)
period <- out$period_seq_all[which.min(out$rss_ls)]
make_plot(lc,period,"panstarrs")
