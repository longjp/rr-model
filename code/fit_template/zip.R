## creates template.zip, an arxiv with all necessary files
## for fitting the model in R and python
rm(list=ls())

direct <- "template"

## remove old zipped info
unlink(direct,recursive=TRUE)
unlink(paste0(direct,".zip"))

## move files
dir.create(direct)
fs <- c("fit_template.R","kappa.py","demo_R.ipynb","demo_python.ipynb",
        "template_sdss.RData","template_des.RData","LC_4099.dat",
        "LC_999886.dat","LC_402316.dat","LC_402316_des.dat","readme")
file.copy(from=fs,to=direct,recursive=TRUE)

## zip and remove folder
zip(paste0(direct,".zip"),direct)
unlink(direct,recursive=TRUE)
