#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from astropy.io import fits
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from astropy.wcs import wcs
import sys
import runsex
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd


# In[ ]:


no_of_scas = 18
no_of_pointings = 44
indexname = '/hpc/group/cosmology/phy-lsst/kxw/akari_long/truth/akari_match_index_Y106_'
fname = '/hpc/group/cosmology/phy-lsst/kxw/akari_long/images/akari_match_Y106_'
snfile = fits.open('/hpc/group/cosmology/phy-lsst/kxw/Matched_Gal_HEAD.FITS')
sndata = snfile[1].data
galfile = fits.open('/hpc/group/cosmology/phy-lsst/kxw/akari_long/truth/akari_match_truth_gal.fits')
galdata = galfile[1].data

sne = dict()
found = np.array([])
notfound = np.array([])
snset = set()
detectset = set()
for sca in range(1, 1+no_of_scas):
    for pointing in range(1,no_of_pointings):
##        file = fname + str(pointing) + '_' + str(sca) + '.fits'
##        detections = runsex.runsex(file)
        file = fname + str(pointing) + '_' + str(sca) + '.delme.txt'
        detections = pd.read_csv(file, skiprows=29, sep='\s+', header=None)
        detectx = np.array(detections[12])
        detecty = np.array(detections[13])
        detectcoords = SkyCoord(detectx, detecty, frame='fk5', unit=u.deg)
        detected = 0
        index = indexname + str(pointing) + '_' + str(sca) + '_sn.fits'
        galindex = indexname + '0_'+ str(sca) + '.fits'
        a = fits.open(index)
        b = fits.open(galindex)
        data = a[1].data
        for i in range(len(data)):
            snid = data[i]['ind']
            thissn = sndata[np.asarray(sndata['snid']).astype(int) == data[i]['ind']]
            if data[i]['mag'] < 30 and thissn['redshift_final'][0] > 0.5:
                gal = b[1].data[b[1].data['ind'] == data[i]['hostid']]
                snx = data[i]['x']
                sny = data[i]['y']
                galxdisp = gal['x'] - snx
                galydisp = gal['y'] - sny
                if len(gal) != 0:
                    sncoords = SkyCoord(data[i]['ra'], data[i]['dec'], frame='fk5', unit=u.rad)
                    galcoords = SkyCoord(gal['ra'], gal['dec'], frame='fk5', unit=u.rad)
                    galsndist = galcoords.separation(sncoords).arcsec
                    if galsndist > 1.5 * galdata[gal['ind'][0]]['size']:
                        if data[i]['mag'] < 24:
                            snset.add(snid)
                        sndists = sncoords.separation(detectcoords).arcsec
                        sndistsort = np.argsort(sndists)
                        if sndists[sndistsort[0]] < 0.16:
                            detected = 1
                        if snid not in sne:
                            thissn = sndata[np.asarray(sndata['snid']).astype(int) == data[i]['ind']]
                            sne[snid] = {'total': 1, 'found': detected, 'peakmag': data['mag'][i]}
                        else:
                            sne[snid]['total'] = sne[snid]['total'] + 1
                            sne[snid]['peakmag'] = np.minimum(sne[snid]['peakmag'],data['mag'][i])
                            sne[snid]['found'] = sne[snid]['found'] + detected
                        if detected == 1:
                            found = np.append(found,data['mag'][i])
                            detectset.add(snid)
                        else:
                            notfound = np.append(notfound,data['mag'][i])
                        detected = 0
        a.close()
        b.close()
        #g = open(fname + str(pointing) + '_' + str(sca) + '.reg', mode='w')
        #for i in range(len(detectx)):
        #    g.write('circle ' + str(detectx[i]) + ' ' + str(detecty[i]) + ' 3.0 # color=red width = 3 \n')
        #g.close()
results = open('SNDetectionResults.txt', mode='a')
results.write(str(sne))
results.close()

plt.hist((found,notfound), histtype = 'step',stacked=True,label=('Images with detections','Images without detections'),bins=30)
plt.xlabel('Supernova magnitude')
plt.ylabel('Number of images')
plt.legend()
plt.savefig('SNFound_HighZ.png')
plt.clf()

bins = np.arange(21,30,0.3)
fbin = np.histogram(found,bins)[0]
nbin = np.histogram(notfound,bins)[0]
tbin = fbin + nbin
ffrac = np.divide(fbin, tbin, out=np.zeros(fbin.shape, dtype=float), where=tbin!=0)
nfrac = np.divide(nbin, tbin, out=np.zeros(nbin.shape, dtype=float), where=tbin!=0)
plt.hist(bins[:-1],bins,weights=ffrac)
plt.ylim(0,1)
plt.xlabel('Supernova magnitude')
plt.ylabel('Detected fraction')
plt.savefig('SNEfficiency_HighZ.png')

print(snset)
print(detectset)
