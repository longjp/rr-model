rm(list=ls())

## see readme for number correspondence

## download 1)
download.file("http://www.astro.washington.edu/users/ivezic/sdss/catalogs/stripe82candidateVar_v1.1.dat.gz",
              "raw/stripe82candidateVar_v1.1.dat.gz",method="wget")
f <- read.table("raw/stripe82candidateVar_v1.1.dat.gz")
cnames <- scan("raw/stripe82candidateVar_v1.1.dat.gz",skip=6,nlines=1,what="char")
colnames(f) <- cnames[2:length(cnames)]
write.table(f,"raw/stripe82candidateVar_v1.1.dat",row.names=FALSE,quote=FALSE)


## download 2)
download.file("http://www.astro.washington.edu/users/ivezic/sdss/catalogs/AllLCs.tar.gz",
              "raw/AllLCs.tar.gz",method="wget")
untar("raw/AllLCs.tar.gz",exdir="raw/AllLCs")


## download 3) and 4)
download.file("http://iopscience.iop.org/0004-637X/708/1/717/suppdata/apj326724t3_mrt.txt",
              "raw/apj326724t3_mrt.txt",method="wget")
download.file("http://iopscience.iop.org/0004-637X/708/1/717/suppdata/apj326724t2_mrt.txt"
              "raw/apj326724t2_mrt.txt",method="wget")






