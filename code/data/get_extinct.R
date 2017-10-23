rm(list=ls())
unlink("radec",recursive=TRUE)
dir.create("radec")
cat <- read.table("raw/stripe82candidateVar_v1.1.dat",header=TRUE)
radec <- cat[,c(2,3)]

nn <- 10000
ii <- 0
while(nn*ii < nrow(cat)){
    fname <- paste0("radec/radec",ii,".txt")
    unlink(fname)
    cat("| ra\t| dec\t|\n",file = fname,append = TRUE)
    cat("| double\t| double\t|\n",file = fname,append = TRUE)
    b <- (nn*ii + 1)
    e <- min((nn*(ii+1)),nrow(cat))
    write.table(radec[b:e,],file=fname,row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    ii <- ii + 1
}



temp <- read.table("raw/apj326724t3_mrt.txt",skip=30) ## get distances for rrlyrae, in different file

extc <- read.table("../make_template/extc.dat")
head(extc)
Rrsdss <- extc[extc[,1]=="SDSS" & extc[,2]=="r",][,3]

temp[1,4]/2.751

ebv*R_r = extinction in r
