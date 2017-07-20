## makes the sdss and des templates
rm(list=ls())

print("making sdss template")
source("make_template_sdss.R")

print("making des template")
source("make_template_des.R")
