rm(list=ls())
load("ps_multi_estimate.RData")

within_x <- function(estimates,truth,thres){
    return(apply(abs((estimates - truth)/truth),1,min) < thres)
}


## lcs with incorrect period estimates
p_est <- p_est_new10
to_use <- !within_x(p_est,periods,0.001)
periods <- periods[to_use]
p_est <- p_est[to_use,]
tms <- tms[to_use]



source("func_saw.R")
ComputeResiduals <- function(m,t,phi,omega){
    X <- cbind(1,Saw(omega*t+phi))
    z <- ComputeBeta(m,t,phi,omega)
    if(z[2]<0){
        z[1] <- mean(m)
        z[2] <- 0
    }
    return(m - X%*%z)
}


N <- 1000
phi_grid <- (1:N)/N

jj <- 1

alg_issue <- rep(0,length(tms))
for(jj in 1:length(tms)){
    print(jj)
    rsst <- rep(0,N)
    rsse <- matrix(0,nrow=N,ncol=5)
    for(ii in 1:N){
        temp <- ComputeResiduals(tms[[jj]]$m,tms[[jj]]$t,phi_grid[ii],1/periods[jj])
        rsst[ii] <- sum(temp^2)
        for(kk in 1:5){
            temp <- ComputeResiduals(tms[[jj]]$m,tms[[jj]]$t,phi_grid[ii],1/p_est[jj,kk])
            rsse[ii,kk] <- sum(temp^2)
        }
    }
    rsse <- max(apply(rsse,2,min))
    alg_issue[jj] <- 1*(min(rsst) < rsse)
}
mean(alg_issue)


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


##pairs(cbind(periods,p_est_sine,p_est_new1,p_est_new5,p_est_new10))
