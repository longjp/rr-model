### compute kappa in fisher von mises dist

import numpy as np
from numpy import random


## arguments
##        phi : vector of phases on [0,1)
## value
##        1-step mle of kappa from
##        fisher von mises dist on unit circle
def kappa(phi):
    x = np.sum(np.cos(2*np.pi*phi))
    y = np.sum(np.sin(2*np.pi*phi))
    rbar = np.sqrt(x*x + y*y) / np.size(phi)
    return (rbar*(2.-rbar**2)) / (1.-rbar**2)

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
