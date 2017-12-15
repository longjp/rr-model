## results go here
rm(list=ls())
source("../../common/funcs.R")
source("../../fit_template/fit_template.R")
source("../funcs.R")
source("../../common/plot_funcs.R")

load("0-fit.RData")
load("../../data/clean/des.RData")

## load templates
load("../../fit_template/template_des.RData")
load("../../fit_template/template_sdss.RData")

## = period estimate analysis
## scatterplot, fraction correct, fraction correct by # epochs

## create dust corrected tms
ebv <- extcr / tem_sdss$dust['r']
tmsc_des <- mapply(DustCorrect,tms_des,ebv,MoreArgs=list(tem=tem_des),SIMPLIFY=FALSE)


period_est <- period_est_des[,1]

## make all plots together
for(ii in 1:length(tmsc_des)){
    tm <- tmsc_des[[ii]]
    lc <- TMtoLC(tm)
    p_est <- period_est[ii]
    omega <- 1 / p_est
    coeffs <- ComputeCoeffs(lc,omega,tem_des,use.dust=FALSE)
    pdf(paste0("figs/",ii,"_one.pdf"),height=8,width=12)
    plotLC(lc,p_est,coeffs,tem_des,main=NULL)
    dev.off()
}

