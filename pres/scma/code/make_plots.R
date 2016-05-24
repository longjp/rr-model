rm(list=ls())



cat <- read.table("stripe82candidateVar_v1.1.dat",header=TRUE)

make_plot <- function(id,period){
    tm <- read.table(paste0("AllLCs/LC_",id,".dat"))

    names(tm) <- c("time","band","mag","sigma")
    bands <- unique(tm$band)
    bands <- as.character(bands)
    bands <- sort(bands)

    ylim <- range(tm$mag)
    xlim <- c(0,1)
    pdf(paste0("../figs/folded_",id,".pdf"),width=12,height=6)
    par(mar=c(5,5,1,1))
    plot(0,0,col=0,ylim=rev(ylim),xlim=xlim,xlab=paste0("Phase (period = ",round(period,2),")"),
         ylab="Magnitude",cex.lab=1.5,xaxs='i')
    for(jj in 1:length(bands)){
        temp <- tm[tm$band==bands[jj],]
        segments((temp$time %% period)/period,temp$mag-temp$sigma,
        (temp$time %% period)/period,temp$mag+temp$sigma,col='grey')
        points((temp$time %% period)/period,temp$mag,col=jj,pch=jj,cex=1.5)
    }
    legend("bottomleft",paste0(bands," Band"),col=1:length(bands),pch=1:length(bands),cex=1.5)
    dev.off()
    
    ylim <- range(tm$mag)
    xlim <- range(tm$time)
    pdf(paste0("../figs/unfolded_",id,".pdf"),width=12,height=6)
    par(mar=c(5,5,1,1))
    plot(0,0,col=0,ylim=rev(ylim),xlim=xlim,xlab="Time",ylab="Magnitude",cex.lab=1.5,xaxs='i')
    for(jj in 1:length(bands)){
        temp <- tm[tm$band==bands[jj],]
        segments(temp$time,temp$mag-temp$sigma,temp$time,temp$mag+temp$sigma,col='grey')
        points(temp$time,temp$mag,col=jj,pch=jj,cex=1.5)
    }
    legend("bottomleft",paste0(bands," Band"),col=1:length(bands),pch=1:length(bands),cex=1.5)
    dev.off()
}


rrlyrae <- read.table("apj326724t2_mrt.txt",skip=42)


id <- 4183016
period <- cat$P[cat$ID==id]
make_plot(id,period)

id <- 4099
period <- rrlyrae[rrlyrae$V1==id,3]
make_plot(id,period)


id <- 7904669
period <- cat$P[cat$ID==id]
make_plot(id,period)

id <- 13350
period <- rrlyrae[rrlyrae$V1==id,3]
make_plot(id,period)
