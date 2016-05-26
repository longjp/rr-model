## parameters for simulation
N <- length(tms) ## number of lcs to analyze
NN <- 10 ## number of newton steps
omegas <- GetFreqs(0.2,1) ## frequency grid
mc.cores <- 12 ## number of processors
topN <- 5 ## number of periods to keep
