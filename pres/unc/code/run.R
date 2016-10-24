rm(list=ls())
unlink("figs",recursive=TRUE)
dir.create("figs")

source("classify.R")
source("training.R")
source("tree-construction.R")
source("cluster.R")

file.copy(from=list.files('figs',full.names=TRUE),to='../pres/figs/',overwrite=TRUE)
