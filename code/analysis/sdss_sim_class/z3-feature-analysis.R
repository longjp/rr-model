rm(list=ls())
load("z2-color_feats.RData")
library(randomForest)

## randomly order rr and nonrr
ix <- sample(1:length(cl))
cl <- cl[ix]
feats <- feats[ix,]
feats_FULL <- feats_FULL[ix,]



## colors for plotting
cols <- c("#00000030","red")
names(cols) <- c("not","rr")



#### MAKE SCATTER PLOT OF ALL colors for well and poorly sampled
## we plot from the .02 to 0.98 quantile of each feature
## where the quantiles are computed from the well sampled data
## for both plots (in order to preserve axis limits)
## WELL SAMPLED
f <- feats_FULL[,!grepl("sd",colnames(feats_FULL))]
ax_range <- apply(f,2,function(x){quantile(x,c(0.02,0.98))})
to_use <- rep(TRUE,nrow(f))
for(ii in 1:ncol(ax_range)){
    to_use <- to_use & (f[,ii] > ax_range[1,ii] & f[,ii] < ax_range[2,ii])
}
to_use[is.na(to_use)] <- FALSE
pairs(f[to_use,],col=cols[cl[to_use]],main="Well Sampled")
dev.new()
## POORLY SAMPLED
f <- feats[,!grepl("sd",colnames(feats))]
to_use <- rep(TRUE,nrow(f))
for(ii in 1:ncol(ax_range)){
    to_use <- to_use & (f[,ii] > ax_range[1,ii] & f[,ii] < ax_range[2,ii])
}
to_use[is.na(to_use)] <- FALSE
pairs(f[to_use,],col=cols[cl[to_use]],main="Poorly Sampled")






## 4 random forest output
## well sample w/ errors
rf.fit <- randomForest(feats_FULL,as.factor(cl))
rf.fit ## 2% error
## well sample w/0 errors
rf.fit <- randomForest(feats_FULL[,2*(0:5)+1],as.factor(cl))
rf.fit ## 3.5% error
## poorly sampled w/ errors
rf.fit <- randomForest(na.roughfix(feats),as.factor(cl))
rf.fit ## 10% error
## poorly sampled w/o errors
rf.fit <- randomForest(na.roughfix(feats)[,2*(0:5)+1],as.factor(cl))
rf.fit ## 18% error


## majority class classifier has error rate
min(table(cl) / length(cl))
