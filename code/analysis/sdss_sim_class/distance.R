rm(list=ls())

## load necessary libraries
library('parallel')
library('multiband')
library('rpart')
library('randomForest')
load("../../fit_template/template.RData")
source("../../fit_template/template.R")
source("../../common/funcs.R")
source("../funcs.R")


## data source
load("../../data/clean/sdss_sim_class.RData")
load("results.RData")
source("../params.R")

fig.dir <- "figs_distance"
unlink(fig.dir,recursive=TRUE)
dir.create(fig.dir)


rrmag <- read.table("rrmag.dat",header=TRUE)
rrmag <- rrmag[rrmag$Sys=="SDSS",]
rrmag <- rrmag[order(rrmag$bnd),]


Nrr <- sum(cl=="rr")
coeffs <- matrix(0,nrow=Nrr,ncol=4)
coeffs_true <- matrix(0,nrow=Nrr,ncol=4)
coeffs_ave <- matrix(0,nrow=Nrr,ncol=4)
lpmed <- log10(median(periods[1:Nrr]))
betas_n_n <- rrmag$c0 + rrmag$p1*(lpmed + 0.2) + rrmag$p2*(lpmed + 0.2)^2
names(betas_n_n) <- names(tem$betas)
for(ii in 1:Nrr){
    print(ii)
    tm <- tms[[ii]]
    lc <- TMtoLC(tm)
    ## using estimated period
    lpmed <- log10(period_est[ii])
    betas_n <- rrmag$c0 + rrmag$p1*(lpmed + 0.2) + rrmag$p2*(lpmed + 0.2)^2
    names(betas_n) <- names(tem$betas)
    tem_temp <- tem
    tem_temp$betas <- betas_n
    omega <- 1/period_est[ii]
    coeffs[ii,] <- ComputeCoeffs(lc,omega,tem_temp)
    ## using true period
    lpmed <- log10(periods[ii])
    betas_n <- rrmag$c0 + rrmag$p1*(lpmed + 0.2) + rrmag$p2*(lpmed + 0.2)^2
    names(betas_n) <- names(tem$betas)
    tem_temp <- tem
    tem_temp$betas <- betas_n
    omega <- 1/periods[ii]
    coeffs_true[ii,] <- ComputeCoeffs(lc,omega,tem_temp)
    ## using estimated periods and fixed betas_n_n
    tem_temp <- tem
    tem_temp$betas <- betas_n_n
    omega <- 1/periods[ii]
    coeffs_ave[ii,] <- ComputeCoeffs(lc,omega,tem_temp)
}





rrlyrae <- read.table("apj326724t3_mrt.txt",skip=30)
names(rrlyrae)[1:5] <- c("ID","ra","dec","ar","d")
rrlyrae$ID <- paste0("LC_",rrlyrae$ID,".dat")

colnames(coeffs) <- c("mu","E[B-V]","a","phi")
rr <- data.frame(names(tms)[1:Nrr],coeffs)
names(rr)[1] <- "ID"

colnames(coeffs_true) <- c("mu_true","E[B-V]_true","a_true","phi_true")
rr_true <- data.frame(names(tms)[1:Nrr],coeffs_true)
names(rr_true)[1] <- "ID"

colnames(coeffs_ave) <- c("mu_ave","E[B-V]_ave","a_ave","phi_ave")
rr_ave <- data.frame(names(tms)[1:Nrr],coeffs_ave)
names(rr_ave)[1] <- "ID"

out <- merge(merge(merge(rrlyrae,rr),rr_true),rr_ave)


plot(log10(out$d),out$mu,ave="Specific")
dev.new()
plot(log10(out$d),out$mu_true,main="True")
dev.new()
plot(log10(out$d),out$mu_ave,main="Ave")
