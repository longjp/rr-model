rm(list=ls())
load("ps_multi_estimate.RData")

print("accuracies for sine")
print(paste("1%:",mean(abs((periods - p_est_sine)/periods) < 0.01)))
print(paste("0.1%:",mean(abs((periods - p_est_sine)/periods) < 0.001)))
print(paste("0.01%:",mean(abs((periods - p_est_sine)/periods) < 0.0001)))

print("accuracies for newton")
print(paste("1%:",mean(abs((periods - p_est_new)/periods) < 0.01)))
print(paste("0.1%:",mean(abs((periods - p_est_new)/periods) < 0.001)))
print(paste("0.01%:",mean(abs((periods - p_est_new)/periods) < 0.0001)))
