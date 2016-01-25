rm(list=ls())
load("estimates.RData")
mean(abs((periods - p_est)/periods) < 0.01)
mean(abs((periods - p_est)/periods) < 0.001)
mean(abs((periods - p_est)/periods) < 0.0001)

