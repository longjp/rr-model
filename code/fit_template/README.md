# RR Lyrae Templates

Fit templates to RR Lyrae variable stars to estimate periods, amplitudes, and distances. The code takes data like this:

![alt text](figs/sdss.png "SDSS Stripe 82 RR Lyrae Light Curve")

and outputs 5 parameters (distance modulus $$\mu$$, dust $E[B-V]$, amplitude, frequency, and phase) which can be used to construct the folded light curve:

![alt text](figs/sdss_fold.png "SDSS Stripe 82 RR Lyrae Folded Light Curve")


The model is designed to work even on very sparsely sampled light curves. For example this Dark Energy Survey RR Lyrae

![alt text](figs/des.png "DES Light Curve")


![alt text](figs/des_fold.png "DES Light Curve")



