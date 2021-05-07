#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from astropy.io import fits
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import pylab
from astropy.wcs import wcs
import sys
import runsex
import astropy.units as u
from astropy.coordinates import SkyCoord


# In[ ]:


no_of_scas = 18
no_of_pointings = 1
indexname = '/hpc/group/cosmology/phy-lsst/kxw/akari_long/truth/akari_match_index_Y106_'
fname = '/hpc/group/cosmology/phy-lsst/kxw/akari_long/images/akari_match_Y106_coadd_six_'
gals = dict()
for sca in range(1, 1+no_of_scas):
    for pointing in range(no_of_pointings):
        file = fname + str(sca) + '_drz.fits'
        detections = runsex.runsex(file)
        detectx = np.array(detections.X_WORLD)
        detecty = np.array(detections.Y_WORLD)
        detectcoords = SkyCoord(detectx, detecty, frame='fk5', unit=u.deg)
        detected = 0
        galindex = indexname + str(pointing) + '_'+ str(sca) + '.fits'
        a = fits.open(galindex)
        data = a[1].data
        galcoords = SkyCoord(data['ra'], data['dec'], frame='fk5', unit=u.rad)
        for i, gal in enumerate(galcoords):
            galid = data[i]['ind']
            if min(gal.separation(detectcoords).arcsec) < 0.3:
                detected = 1
            if galid not in gals:
                gals[galid] = {'total': 1, 'found': detected, 'galmag': data[i]['mag']}
            else:
                gals[galid]['total'] = gals[galid]['total'] + 1
                gals[galid]['found'] = gals[galid]['found'] + detected
            detected = 0
        a.close()
results = open('CoaddResults.txt', mode='w')
results.write(str(gals))
results.close()
found = np.array([])
notfound = np.array([])
for v in gals.values():
    if float(v['found']) / float(v['total']) >= 0.5:
        found = np.append(found, v['galmag'])
    else:
        notfound = np.append(notfound, v['galmag'])
bins = np.arange(17,28,1)
fbin = np.histogram(found,bins)[0]
nbin = np.histogram(notfound,bins)[0]
tbin = fbin + nbin
ffrac = fbin/tbin
nfrac = nbin/tbin
plt.hist(bins[:-1],bins,weights=ffrac)
plt.ylim(0,1)
plt.savefig('CoaddFound.png')
plt.clf()
plt.hist(bins[:-1],bins,weights=nfrac)
plt.ylim(0,1)
plt.savefig('CoaddNotFound.png')
plt.clf()
plt.hist(found,range=(15,30),histtype='step',bins=30)
plt.hist(np.append(found,notfound),range=(15,30),histtype='step',bins=30)
plt.savefig('CoaddComp.png')
