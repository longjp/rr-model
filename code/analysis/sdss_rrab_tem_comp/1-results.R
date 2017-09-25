rm(list=ls())

## load necessary libraries
library('parallel')
load("../../fit_template/template_sdss.RData")
source("../../fit_template/fit_template.R")
source("../../common/funcs.R")
source("../funcs.R")

## data source
load("../../data/clean/sdss_sim_class.RData")

## parameters for simulation
##source("../params.R")

## load results from sim
load("0-fit.RData")

## relies on fact that RRL are first in list
N <- sum(cl=="rr")
N <- nrow(period_est_new)


periods <- periods[1:N]


## make 3 pairs of scatterplots of
## 1. true periods vs. old/new
## 2. sesar distances vs. old/new
## 3. schlegel dust vs. old/new

