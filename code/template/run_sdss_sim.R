## run simulation
rm(list=ls())
unlink("*.RData")

rm(list=ls())
source("make_template.R")

rm(list=ls())
source("sdss_sim.R")

rm(list=ls())
source("estimate_params.R")

rm(list=ls())
source("estimate_params_lomb.R")

rm(list=ls())
source("results.R")
