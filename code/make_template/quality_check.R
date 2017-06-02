### fits well sampled light curves to templates, calculates median residuals
### used to check that templates are working well, compare minor
### implementation differences
rm(list=ls())
source('../common/funcs.R')
source('../fit_template/fit_template.R')
load("../data/clean/sdss_rrab.RData")
load("../fit_template/template.RData")


## find model error
med_res <- vector("numeric",length(tms))
for(ii in 1:length(med_res)){
    lc <- TMtoLC(tms[[ii]])
    coeffs <- ComputeCoeffs(lc,1/periods[ii],tem)
    preds <- PredictTimeBand(lc[,1],lc[,2],1/periods[ii],coeffs,tem)
    med_res[ii] <- median(abs(preds - lc[,3]))
}


print("model error:")
hist(med_res)
median(med_res) ## .031 seems reasonable for this number


## find model error caused by shape
med_res <- vector("numeric",length(tms))
for(ii in 1:length(med_res)){
    temp <- tms[[ii]]
    for(jj in 1:length(temp)){
        temp[[jj]][,2] <- temp[[jj]][,2] - mean(temp[[jj]][,2]) + tem$betas[jj]
    }
    lc <- TMtoLC(temp)
    coeffs <- ComputeCoeffs(lc,1/periods[ii],tem)
    preds <- PredictTimeBand(lc[,1],lc[,2],1/periods[ii],coeffs,tem)
    med_res[ii] <- median(abs(preds - lc[,3]))
}

print("model error, perfect mean:")
hist(med_res)
median(med_res) ## perfect mean offers only small advantage over actual model


## how much error due only to photometry?
med_res <- vector("numeric",length(tms))
for(ii in 1:length(med_res)){
    lc <- TMtoLC(tms[[ii]])
    med_res[ii] <- median(abs(rnorm(nrow(lc),mean=0,sd=lc[,4])))
}

print("photometric error:")
hist(med_res)
median(med_res) ## perfect mean offers only small advantage over actual model


## are templates mean 0?
print("the mean of each template is:")
rowMeans(tem$templates)




## errors are
lcs <- lapply(tms,TMtoLC)
sds <- unlist(lapply(lcs,function(x){x[,4]}))
summary(sds*.6)
