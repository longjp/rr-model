## moves files into separate template fitting only repo
## for uploading to github
rm(list=ls())


args = commandArgs(trailingOnly=TRUE)

direct <- args[1]
print(paste0("copying files to: ",direct))

## dir.create(direct,showWarnings=FALSE)
## fs <- list.files(direct,full.names=TRUE)
## unlink(fs,recursive=TRUE)

## move files
fs <- c("fit_template.R","demo_R.ipynb","demo_python.ipynb",
        "template_sdss.RData","template_des.RData","LC_4099.dat",
        "LC_999886.dat","LC_402316.dat","LC_402316_des.dat","README.md",
        "figs/","LICENSE")
print(fs)
print(direct)
file.copy(from=fs,to=direct,recursive=TRUE)

## zip and remove folder
##zip(paste0(direct,".zip"),direct)
##unlink(direct,recursive=TRUE)
