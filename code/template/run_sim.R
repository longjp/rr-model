## run simulation
rm(list=ls())
unlink("*.RData")

rm(list=ls())
source("make_template.R")

rm(list=ls())
source("sim.R")

rm(list=ls())
source("estimate_params.R")

rm(list=ls())
source("results.R")
