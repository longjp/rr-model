rm(list=ls())
unlink("radec",recursive=TRUE)
dir.create("radec")
cat <- read.table("../raw/stripe82candidateVar_v1.1.dat",header=TRUE)
radec <- cat[,c(2,3)]

nn <- 10000
ii <- 0
while(nn*ii < nrow(cat)){
    fname <- paste0("radec/radec",ii,".dat")
    unlink(fname)
    cat("| ra\t| dec\t|\n",file = fname,append = TRUE)
    cat("| double\t| double\t|\n",file = fname,append = TRUE)
    b <- (nn*ii + 1)
    e <- min((nn*(ii+1)),nrow(cat))
    write.table(radec[b:e,],file=fname,row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    ii <- ii + 1
}
