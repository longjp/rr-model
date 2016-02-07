rm(list=ls())
load("ps_multi_estimate.RData")

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


##pairs(cbind(periods,p_est_sine,p_est_new1,p_est_new5,p_est_new10))
