## load des RRL lightcurves in Stripe 82 region
## associated light curves with DES filename and sloan ID
rm(list=ls())
source('../common/funcs.R')

## read in katelyn catalog data and get sdss ids in nice form
cat <- read.table("raw/known_rr_lcs/rrab_only.tab",
                  colClasses=c("numeric","character",
                               "character","character","numeric"),
                  header=TRUE)
cat$ID <- as.character(gsub(".0","",cat$ID,fixed=TRUE))
cat$filename <- gsub("./","",cat$filename,fixed=TRUE)
sdss_id <- cat$ID

## read in des lcs
lcs_des <- vector("list",nrow(cat))
for(ii in 1:length(lcs_des)){
    fname <- gsub("./","",cat$filename[ii],fixed=TRUE)
    lc <- read.table(paste0("raw/known_rr_lcs/",fname),header=TRUE,stringsAsFactors=FALSE)
    lc <- lc[,c(1,4,2,3)]
    names(lc) <- c("time","band","mag","error")
    lcs_des[[ii]] <- lc
}
tms_des <- lapply(lcs_des,LCtoTM)
names(tms_des) <- cat$filename

## save output
save(tms_des,sdss_id,file="clean/des.RData")
