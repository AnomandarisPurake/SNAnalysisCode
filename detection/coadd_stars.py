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
import glob
import pandas as pd

# In[ ]:


no_of_scas = 18
no_of_pointings = 1
indexname = '../truth/paper_final_index_Y106_'
gals = dict()

def truth_gals(ra_min, ra_max, dec_min, dec_max):
    garbage = fits.open(indexname + '132_1.fits')
    garbage_length = len(garbage[1].data)
    gal_list = garbage[1].data
    for sca in range(1, 1+no_of_scas):
        for pointing in range(132,176):
            galindex = indexname + str(pointing) + '_'+ str(sca) + '.fits'
            a = fits.open(galindex)
            data = a[1].data
            true_gals = data[(data['ra'] > ra_min) & (data['ra'] < ra_max) & (data['dec'] > dec_min) & (data['dec'] < dec_max)]
            gal_list = np.append(gal_list, true_gals)
            a.close()
    return gal_list[garbage_length:]
    
def coadd_truth(filename):
    name = '../truth/paper_final_index_' + filename[-15:-4] + '.fits.gz'
    gals = fits.open(name)
    data = gals[1].data
    data = data[data['gal_star'] == 1]
    #print(len(data[data['mag'][:,3] != 0]))
    return data[data['mag'][:,3] != 1]


for f in glob.glob('paper_final_Y*.cat'):
    detected = 0
    cfile = f
    cdetections = pd.read_csv(cfile, skiprows=54, sep='\s+', header=None)
    cdetectx = np.array(cdetections[46])
    cdetecty = np.array(cdetections[47])
    cdetectra = np.array(cdetections[48])
    cdetectdec = np.array(cdetections[49])
    cdetectcoords = SkyCoord(cdetectra, cdetectdec, frame='fk5', unit=u.deg)
    ra_min = min(cdetectra) * np.pi / 180
    ra_max = max(cdetectra) * np.pi / 180
    dec_min = min(cdetectdec) * np.pi / 180
    dec_max = max(cdetectdec) * np.pi / 180
    #bdata = truth_gals(ra_min, ra_max, dec_min, dec_max)
    #print(len(bdata))
    data = coadd_truth(f)
    #print(len([x for x in bdata if x not in data['ind']]))
    #print(len([x for x in data['ind'] if x not in bdata])) 
    #for x in bdata:
        #match = data[data['ind'] == x['ind']]
        #if match:
            #print('match')
            #print(match)
            #print('single_truth')
            #print(x['ra'] * 180 / np.pi, x['dec'] * 180 / np.pi)
            #print(x['mag'])
            #print('coadd_truth')
            #print(match['mag'])
            #print(match['ra'], match['dec'])
    galcoords = SkyCoord(data['ra'], data['dec'], frame='fk5', unit=u.deg)
    for i, gal in enumerate(galcoords):
        galid = data[i]['ind']
        if min(gal.separation(cdetectcoords).arcsec) < 1:
            detected = 1
        if galid not in gals:
            gals[galid] = {'total': 1, 'found': detected, 'galmag': data[i]['mag'][3]}
        else:
            gals[galid]['total'] = gals[galid]['total'] + 1
            gals[galid]['found'] = gals[galid]['found'] + detected
        detected = 0

results = open('CoaddStarResults.txt', mode='w')
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
plt.ylabel('Detection fraction')
plt.xlabel('Y Band Magnitude')
plt.savefig('CoaddStarFound.png')
plt.clf()
