## python wrapper for template.R
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects.vectors import IntVector, FloatVector, StrVector
import csv
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

r=robjects.r
r.source("template.R")
r.load("template.RData")
FitTemplate=robjects.r['FitTemplate']

## load a light curve
fname="LC_4099.dat"
with open(fname) as csvf:
    f = csv.reader(csvf,delimiter=' ')
    time, band, mag, error = zip(*f)

time = np.array(time,dtype='float64')
mag = np.array(mag,dtype='float64')
error = np.array(error,dtype='float64')

## create R dataframe using time,band,mag,error
lc = robjects.r['TBMEtoLC'](FloatVector(time),StrVector(band),FloatVector(mag),FloatVector(error))

## choose frequency grid
omegas=FloatVector(np.arange(start=1.0,stop=5.0,step=0.1/4000.0))

## compute rss (residual sum of squares) for each frequency
rss=FitTemplate(lc,omegas,r.tem)
## select best fitting period
pest=1.0/omegas[np.argmin(rss)]

## plot lightcurve folded on estimated period
cols = {'g': 'green', 'i': 'red', 'r':'blue','u': 'orange', 'z': 'yellow'}
pts=plt.scatter(np.mod(time,pest),mag,color=map(cols.get,band))
plt.gca().invert_yaxis()
plt.xlim([0,pest])
class_colours = map(cols.get,band)
recs = []
for i in range(0,len(class_colours)):
    recs.append(mpatches.Rectangle((0,0),1,1,fc=class_colours[i]))

plt.legend(recs,cols.keys(),loc=4)
plt.show()



