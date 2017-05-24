### fits well sampled light curves to templates, calculates median residuals
### used to check that templates are working well, compare minor
### implementation differences
rm(list=ls())
source('../common/funcs.R')
source('../fit_template/fit_template.R')
load("../data/clean/sdss_rrab.RData")
load("../fit_template/template.RData")

med_res <- vector("numeric",length(tms))
for(ii in 1:length(med_res)){
    lc <- TMtoLC(tms[[ii]])
    coeffs <- ComputeCoeffs(lc,1/periods[ii],tem)
    preds <- PredictTimeBand(lc[,1],lc[,2],1/periods[ii],coeffs,tem)
    med_res[ii] <- median(abs(preds - lc[,3]))
}


hist(med_res)
median(med_res) ## .031 seems reasonable for this number


rowMeans(tem$templates)
