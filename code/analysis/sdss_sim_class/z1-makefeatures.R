## compute uncertainties on all parameter estimates
## output sds (square root diagonal of hessian) along with features

require(numDeriv)
load("../../fit_template/template_sdss.RData")
load("../../fit_template/template_des.RData")
source("../../fit_template/fit_template.R")


ComputeUncertainty <- function(coeffs,omega,lc,tem,use.errors=TRUE,use.dust=TRUE){
    if(use.dust){
        use.dust <- CheckNumberBands(lc)
    }
    CheckLC(lc)
    ## unpack coefficients
    tem <- CheckTemLC(tem,lc)
    dat <- AugmentData(lc,tem,use.errors)
    mag <- dat[[1]]$mag
    dust <- dat[[1]]$dust
    t <- dat[[1]]$time
    weights <- 1 / dat[[1]]$error^2
    nb <- dat[[2]]
    coeffs <- c(omega,coeffs)
    return(solve(hessian(ComputeRSS,coeffs,mag=mag,tem=tem,nb=nb,t=t,dust=dust,weights=weights)))
}

ComputeRSS <- function(coeffs,mag,tem,nb,t,dust,weights){
    omega <- coeffs[1]
    phi <- coeffs[5]
    coeffs <- coeffs[2:4]
    mag <- mag - rep.int(tem$abs_mag(1/omega,tem)[1,],nb)
    gammaf <- ConstructGamma(t,nb,phi,omega,tem$template_funcs)
    resid <- mag - coeffs[1] - coeffs[2]*dust - coeffs[3]*gammaf
    return(sum(weights*resid^2))
}





## read in des data for same light curve
fname <- "LC_402316_des.dat"
lc_des <- read.table(fname,header=TRUE,stringsAsFactors=FALSE)
lc_des <- lc_des[,c(1,4,2,3)]
names(lc_des) <- c("time","band","mag","error")

## construct frequency grid
freq_min <- 1.0/0.95 ## max period of rrab is about 0.95 days 
freq_max <- 1.0/0.4 ## min period of rrab is about 0.4 days
max_phase_error <- .01 ## maximum phase error fraction, see documentation
freq_space <- (max_phase_error*4)/(max(lc_des[,1]) - min(lc_des[,1]))
omegas <- seq(from=freq_min,to=freq_max,by=freq_space)


## estimate period with des templates
rss_des <- FitTemplate(lc_des,omegas,tem_des)
omega_est_des <- omegas[which.min(rss_des)]
p_est_des <- 1/omega_est_des
coeffs_des <- ComputeCoeffs(lc_des,omega_est_des,tem_des)
names(coeffs_des) <- c("mu","ebv","amp","phase")


des_poor <- ComputeUncertainty(coeffs_des,omega_est_des,lc_des,tem_des)
des_poor
