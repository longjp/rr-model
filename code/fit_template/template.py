## python wrapper for template.R
import numpy as np
import rpy2.robjects as robjects
import csv
import matplotlib.pyplot as plt
from rpy2.robjects.vectors import IntVector, FloatVector, StrVector
r=robjects.r
r.source("template.R")
r.load("template.RData")

## load a light curve
fname="LC_4099.dat"
with open(fname) as csvf:
    f = csv.reader(csvf,delimiter=' ')
    time, band, mag, error = zip(*f)

time = np.array(time,dtype='float64')
mag = np.array(mag,dtype='float64')
error = np.array(error,dtype='float64')

## create R lc dataframe using time,band,mag,error
lc = robjects.r['TBMEtoLC'](FloatVector(time),StrVector(band),FloatVector(mag),FloatVector(error))

## choose frequency grid
omegas=FloatVector(np.arange(start=1.0,stop=5.0,step=0.1/4000.0))

## find rss
FitTemplate=robjects.r['FitTemplate']
rss=FitTemplate(lc,omegas,r.tem)
pest=1.0/omegas[np.argmin(rss)]

## plot lightcurve and folded
plt.scatter(np.mod(time,pest),mag)
plt.show()


## TODO:
## invert y axis
## color points by bands
## plot error bars
## is R code using errors
