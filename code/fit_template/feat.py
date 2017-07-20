### compute kappa in fisher von mises dist

import numpy as np
from numpy import random

def kappa(phi):
    n = np.size(phi)
    sumcos = np.sum(np.cos(2*np.pi*phi))
    sumsin = np.sum(np.sin(2*np.pi*phi))
    rbar = np.sqrt(sumcos*sumcos + sumsin*sumsin) / n
    rbar2 = np.power(rbar,2)
    return (rbar*(2.-rbar2)) / (1.-rbar2)


## generate some random data on [0,1)
mu = 0.5
sigma = .1
phi = np.random.normal(mu, sigma, 1000)
phi = phi % 1.0
kappa(phi)


## generate some random data on [0,1)
mu = 0.9
sigma = .1
phi = np.random.normal(mu, sigma, 1000)
phi = phi % 1.0
kappa(phi)


## generate some random data on [0,1)
mu = 0.5
sigma = .02
phi = np.random.normal(mu, sigma, 1000)
phi = phi % 1.0
kappa(phi)
     
