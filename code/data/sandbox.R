rm(list=ls())


rrlyrae <- read.table("raw/apj326724t3_mrt.txt",skip=30)


names(rrlyrae)[1:5] <- c("ID","ra","dec","ar","d")

plot(rrlyrae$ra,rrlyrae$dec)
plot(rrlyrae$ra,rrlyrae$dec)


plot(rrlyrae$d*cos(2*pi*(rrlyrae$ra)/360),rrlyrae$d*sin(2*pi*(rrlyrae$ra)/360))







rrlyrae <- read.table("raw/apj326724t2_mrt.txt",skip=42)

cat <- read.table("raw/stripe82candidateVar_v1.1.dat",header=TRUE)
names(rrlyrae)[1] <- "ID"

rrlyrae <- merge(cat,rrlyrae)

plot(rrlyrae$ra,rrlyrae$dec)
