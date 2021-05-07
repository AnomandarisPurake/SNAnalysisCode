#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import wcs
import sys
import runsex
import plotsetup


# In[ ]:


#Single Image General Galaxy Results
plotsetup.fullpaperfig()
gals = eval(open('GalResults.txt',mode='r').read())
found = np.array([])
notfound = np.array([])
for v in gals.values():
    if float(v['found']) / float(v['total']) >= 0.5:
        found = np.append(found, v['galmag'])
    else:
        notfound = np.append(notfound, v['galmag'])
bins = np.arange(17,30,1)
fbin = np.histogram(found,bins)[0]
nbin = np.histogram(notfound,bins)[0]
tbin = fbin + nbin
ffrac = fbin/tbin
nfrac = nbin/tbin
plt.hist(bins[:-1],bins,weights=ffrac,histtype='step',label='Single Images')
plt.ylim(0,1)

#Coadd General Galaxy Results
gals = eval(open('CoaddResults.txt',mode='r').read())
found = np.array([])
notfound = np.array([])
for v in gals.values():
    if float(v['found']) / float(v['total']) >= 0.5:
        found = np.append(found, v['galmag'])
    else:
        notfound = np.append(notfound, v['galmag'])
bins = np.arange(17,30,1)
fbin = np.histogram(found,bins)[0]
nbin = np.histogram(notfound,bins)[0]
tbin = fbin + nbin
ffrac = fbin/tbin
nfrac = nbin/tbin
plt.hist(bins[:-1],bins,weights=ffrac,histtype='step',label='Coadded Images')
plt.legend()
plt.savefig('GalaxyEfficiency.png')
plt.clf()

    

