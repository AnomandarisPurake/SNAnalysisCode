#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import wcs
import sys
import runsex


# In[ ]:


no_of_scas = 18
no_of_pointings = 44
indexname = '/hpc/group/cosmology/phy-lsst/kxw/akari_long/truth/akari_match_index_Y106_'
fname = '/hpc/group/cosmology/phy-lsst/kxw/akari_long/images/akari_match_Y106_'
gals = dict()
for sca in range(1, 1+no_of_scas):
    for pointing in range(no_of_pointings):
        file = fname + str(pointing) + '_' + str(sca) + '.fits'
        detections = runsex.runsex(file)
        detectx = detections.X_IMAGE
        detecty = detections.Y_IMAGE
        detected = 0
        galindex = indexname + str(pointing) + '_'+ str(sca) + '.fits'
        a = fits.open(galindex)
        data = a[1].data
        for i in range(len(data)):
            galx = data[i]['x']
            galy = data[i]['y']
            galid = data[i]['ind']
            for j in range(len(detectx)):
                if (detectx[j] - galx)**2 + (detecty[j] - galy)**2 < 3:
                    detected = 1
            if galid not in gals:
                gals[galid] = {'total': 1, 'found': detected, 'galmag': data[i]['mag']}
            else:
                gals[galid]['total'] = gals[galid]['total'] + 1
                gals[galid]['found'] = gals[galid]['found'] + detected
            detected = 0
        a.close()
results = open('GalResults.txt', mode='w')
results.write(str(gals))
results.close()
found = np.array([])
notfound = np.array([])
for v in gals.values():
    if float(v['found']) / float(v['total']) >= 0.5:
        found = np.append(found, v['galmag'])
    else:
        notfound = np.append(notfound, v['galmag'])
bins = np.arange(17,26,1)
fbin = np.histogram(found,np.arange(17,26,1))[0]
nbin = np.histogram(notfound,np.arange(17,26,1))[0]
tbin = fbin + nbin
ffrac = fbin/tbin
nfrac = nbin/tbin
plt.hist(bins[:-1],bins,weights=ffrac)
plt.ylim(0,1)
plt.savefig('GalFound.png')
plt.clf()
plt.hist(bins[:-1],bins,weights=nfrac)
plt.ylim(0,1)
plt.savefig('GalNotFound.png')

    

