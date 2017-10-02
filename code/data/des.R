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


## use des light curves with sdss crossmatch
## want to order light curves by filename
load("clean/sdss_rrab.RData")
tms_sdss <- tms
names(tms_sdss) <- gsub(".dat","",gsub("LC_","",names(tms_sdss)))
to_use <- sdss_id %in% names(tms_sdss)
sum(to_use)
tms_des <- tms_des[to_use]
sdss_id <- sdss_id[to_use]
tms_des <- tms_des[order(sdss_id)]

## get sdss sources which are in des, order
to_use <- names(tms_sdss) %in% sdss_id
sum(to_use)
tms_sdss <- tms_sdss[to_use]
periods <- periods[to_use]
distance <- distance[to_use]
ra <- ra[to_use]
dec <- dec[to_use]
extcr <- extcr[to_use]

ords <- order(names(tms_sdss))
tms_sdss <- tms_sdss[ords]
periods <- periods[ords]
distance <- distance[ords]
ra <- ra[ords]
dec <- dec[ords]
extcr <- extcr[ords]

## save output
save(tms_des,tms_sdss,periods,distance,ra,dec,extcr,file="clean/des.RData")
