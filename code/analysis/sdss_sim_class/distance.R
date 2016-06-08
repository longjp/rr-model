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


period_est <- period_est[,1]

Nrr <- sum(cl=="rr")
coeffs_model <- matrix(0,nrow=Nrr,ncol=4)
coeffs_lomb <- matrix(0,nrow=Nrr,ncol=4)
for(ii in 1:Nrr){
    tm <- tms[[ii]]
    lc <- TMtoLC(tm)
    ## distance estimates with rr model
    omega <- 1/period_est[ii]
    coeffs_model[ii,] <- ComputeCoeffs(lc,omega,tem)
    ## distance estimates using lomb
    omega <- 1/period_est_lomb[ii]
    coeffs_lomb[ii,] <- ComputeCoeffs(lc,omega,tem)
}


rrlyrae <- read.table("apj326724t3_mrt.txt",skip=30)
names(rrlyrae)[1:5] <- c("ID","ra","dec","ar","d")
rrlyrae$ID <- paste0("LC_",rrlyrae$ID,".dat")

colnames(coeffs_model) <- c("mu_model","E[B-V]_model","a_model","phi_model")
rr_model <- data.frame(names(tms)[1:Nrr],coeffs_model)
names(rr_model)[1] <- "ID"

colnames(coeffs_lomb) <- c("mu_lomb","E[B-V]_lomb","a_lomb","phi_lomb")
rr_lomb <- data.frame(names(tms)[1:Nrr],coeffs_lomb)
names(rr_lomb)[1] <- "ID"


out <- Reduce(function(x, y) merge(x, y), list(rrlyrae,rr_model,rr_lomb))

pdf("distance_comparison.pdf")
par(mar=c(5,5,1,1))
plot(out$d,10^(out$mu_model/5 + 1)/1000,
     xlab="Ground Truth Distance in kpc (Sesar 2010)",
     ylab="Estimate in kpc",
     cex.lab=1.3)
abline(a=0,b=1)
dev.off()
