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
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd

# In[ ]:


no_of_scas = 18
no_of_pointings = 1
indexname = '../truth/paper_final_index_Y106_'
fname = 'paper_final_Y106_'
galfile = fits.open('../truth/paper_final_truth_gal.fits')
galdata = galfile[1].data
gals = dict()
for sca in range(1, 1+no_of_scas):
    for pointing in range(132, 176):
        file = fname + str(pointing) + '_' + str(sca) + '.cat'
        detections = pd.read_csv(file, skiprows=54, sep='\s+', header=None)
        detectx = np.array(detections[48])
        detecty = np.array(detections[49])
        detectcoords = SkyCoord(detectx, detecty, frame='fk5', unit=u.deg)
        detected = 0
        galindex = indexname + str(pointing) + '_'+ str(sca) + '.fits'
        a = fits.open(galindex)
        data = a[1].data
        galcoords = SkyCoord(data['ra'], data['dec'], frame='fk5', unit=u.rad)
        for i, gal in enumerate(galcoords):
            galid = data[i]['ind']
#            galaxy = galdata[galdata['gind'] == galid]
            if min(gal.separation(detectcoords).arcsec) < 1:
                detected = 1
            if galid not in gals:
                gals[galid] = {'total': 1, 'found': detected, 'galmag': data[i]['mag']}
            else:
                gals[galid]['total'] = gals[galid]['total'] + 1
                gals[galid]['found'] = gals[galid]['found'] + detected
            detected = 0
        a.close()

results = open('SingleResults.txt', mode='w')
results.write(str(gals))
results.close()

found = np.array([])
notfound = np.array([])
for v in gals.values():
    if float(v['found']) / float(v['total']) >= 0.5:
        found = np.append(found, v['galmag'])
    else:
        notfound = np.append(notfound, v['galmag'])
bins = np.arange(16,28,0.5)
fbin = np.histogram(found,bins)[0]
nbin = np.histogram(notfound,bins)[0]
tbin = fbin + nbin
ffrac = fbin/tbin
nfrac = nbin/tbin
plt.hist(bins[:-1],bins,weights=ffrac,histtype='step')
plt.ylim(0,1)
plt.xlabel('Y Band Mag')
plt.ylabel('Detected fraction')
plt.savefig('SingleFound.png')
plt.clf()
