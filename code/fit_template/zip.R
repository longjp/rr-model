## creates template.zip, an arxiv with all necessary files
## for fitting the model in R and python
rm(list=ls())

direct <- "template"

## remove old zipped info
unlink(direct,recursive=TRUE)
unlink(paste0(direct,".zip"))

## move files
dir.create(direct)
fs <- c("fit_template.R","feat.py","template.R","template.py",
        "template_sdss.RData","template_des.RData","LC_4099.dat",
        "LC_999886.dat","readme_zip")
file.copy(from=fs,to=direct,recursive=TRUE)

file.rename(from=paste0(direct,"/readme_zip"),
            to=paste0(direct,"/readme"))

## zip and remove folder
zip(paste0(direct,".zip"),direct)
unlink(direct,recursive=TRUE)
