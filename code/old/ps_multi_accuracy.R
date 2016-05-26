rm(list=ls())
source('func_multi_saw.R')
load("ps_multi_estimate.RData")
load("ps_multi.RData")


within_x <- function(estimates,truth,thres){
    return(apply(abs((estimates - truth)/truth),1,min) < thres)
}

print("accuracies for sine")
print(paste("1%:",mean(within_x(p_est_sine,periods,0.01))))
print(paste("0.1%:",mean(within_x(p_est_sine,periods,0.001))))
print(paste("0.01%:",mean(within_x(p_est_sine,periods,0.0001))))
print("")

print("accuracies for newton: NN=1")
print(paste("1%:",mean(within_x(p_est_new1,periods,0.01))))
print(paste("0.1%:",mean(within_x(p_est_new1,periods,0.001))))
print(paste("0.01%:",mean(within_x(p_est_new1,periods,0.0001))))
print("")

print("accuracies for newton: NN=5")
print(paste("1%:",mean(within_x(p_est_new5,periods,0.01))))
print(paste("0.1%:",mean(within_x(p_est_new5,periods,0.001))))
print(paste("0.01%:",mean(within_x(p_est_new5,periods,0.0001))))
print("")

print("accuracies for newton: NN=10")
print(paste("1%:",mean(within_x(p_est_new10,periods,0.01))))
print(paste("0.1%:",mean(within_x(p_est_new10,periods,0.001))))
print(paste("0.01%:",mean(within_x(p_est_new10,periods,0.0001))))
print("")

p_est <- p_est_new10[,1]


for(ii in 1:length(p_est)){
    lc <- tms[[ii]]
    lc_band <- tms_band[[ii]]
    omega <- 1/p_est[ii]
    params <- GetParamsPhiGrid(lc$m,lc$t,lc$ampj,lc$phij,omega,cc)
    pdf(paste0("figs/",ii,".pdf"),height=12,width=8)
    band_names <- unique(lc_band)
    N <- 1000
    par(mar=c(5,3,3,1),mfcol=c(5,1))
    for(jj in 1:length(band_names)){
        lc_temp <- lc[lc_band==band_names[jj],]
        t <- ((1:N)/N)/omega
        ampj <- lc_temp$ampj[1]
        phij <- lc_temp$phij[1]
        pred <- params["beta"] + params["amp"]*ampj*Saw(omega*(t+phij)+params["phi"],cc)
        plot(t * omega,pred,type='l',xlab="Phase",ylab="Mag",main=band_names[jj])
        points((lc_temp$t %% (1/omega))*omega,lc_temp$m)
    }
    dev.off()
}

##pairs(cbind(periods,p_est_sine,p_est_new1,p_est_new5,p_est_new10))
