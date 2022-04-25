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
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd
from scipy import stats


# In[ ]:


no_of_scas = 18
no_of_pointings = 44
indexname = '../truth/paper_final_index_Y106_'
fname = 'paper_final_Y106_'
stars = dict()
found = np.array([])
notfound = np.array([])

for sca in range(1, 1+no_of_scas):
    for pointing in range(132,176):
##        file = fname + str(pointing) + '_' + str(sca) + '.fits'
##        detections = runsex.runsex(file)
        file = fname + str(pointing) + '_' + str(sca) + '.cat'
        detections = pd.read_csv(file, skiprows=54, sep='\s+', header=None)
        detectx = np.array(detections[48])
        detecty = np.array(detections[49])
        detectcoords = SkyCoord(detectx, detecty, frame='fk5', unit=u.deg)
        detectmag = np.array(detections[15])
        detectmagerr = np.array(detections[18])
        detected = 0
        index = indexname + str(pointing) + '_' + str(sca) + '_star.fits'
        a = fits.open(index)
        data = a[1].data
        for i in range(len(data)):
            starid = data[i]['ind']
            if data[i]['mag'] < 30:
                starx = data[i]['x']
                stary = data[i]['y']
                starcoords = SkyCoord(data[i]['ra'], data[i]['dec'], frame='fk5', unit=u.rad)
                stardists = starcoords.separation(detectcoords).arcsec
                stardistsort = np.argsort(stardists)
                if stardists[stardistsort[0]] < 0.16:
                    detected = 1
                if starid not in stars:
                    stars[starid] = {'total': 1, 'found': detected, 'mag': data[i]['mag'], 'Smag': detectmag[stardistsort[0]], 'SNR': 1.0/detectmagerr[stardistsort[0]]}
                else:
                    stars[starid]['total'] = stars[starid]['total'] + 1
                    if detected == 1:
                        stars[starid]['SNR'] = stars[starid]['SNR'] + 1.0/detectmagerr[stardistsort[0]]
                    stars[starid]['found'] = stars[starid]['found'] + detected
                detected = 0
        a.close()

results = open('StarSNRResults.txt', mode='w')
results.write(str(stars))
results.close()

SNR = np.array([])
mag = np.array([])
Smag = np.array([])
for v in stars.values():
    if v['found'] >= 10:
        SNR = np.append(SNR, v['SNR']/v['found'])
        mag = np.append(mag, v['mag'])
        Smag = np.append(Smag, v['Smag'])

plt.scatter(mag,Smag,s = 4)
plt.xlabel('Truth mag')
plt.ylabel('Sextractor mag')
plt.savefig('MagComparison.png')
plt.clf()

plt.scatter(mag,SNR,s=4,alpha=0.1)
plt.xlabel('Star magnitude')
plt.ylabel('SNR')
plt.savefig('StarSNR.png')
plt.clf()

bin_mean, bin_mean_edge, bin_mean_number = stats.binned_statistic(mag,SNR,bins=20,range=(20,30))
bin_median, bin_median_edge, bin_median_number = stats.binned_statistic(mag,SNR,statistic='median',bins=20,range=(20,30))
plt.hlines(bin_mean, bin_mean_edge[:-1], bin_mean_edge[1:])
plt.xlabel('Binned star magnitude')
plt.ylabel('Mean SNR')
plt.savefig('MeanStarSNR.png')
plt.clf()

plt.hlines(bin_median, bin_median_edge[:-1], bin_median_edge[1:])
plt.xlabel('Binned star magnitude')
plt.ylabel('Median SNR')
plt.savefig('MedianStarSNR.png')
plt.clf()
    
other = pd.read_csv('OUT_FOR_KEVIN.OUTLIER.TEXT',skiprows=15, sep='\s+',header=None)
logSNR = np.array(other[17])
fluxcal = np.array(other[15])
band = np.array(other[5])
field = np.array(other[3])
deep = np.intersect1d(np.asarray(band == 'Y').nonzero()[0],np.asarray(field == 'DEEP').nonzero()[0])
deep_mag = 27.5 - 2.5 * np.log10(fluxcal[deep])
deep_SNR = np.power(10,logSNR[deep])
plt.scatter(mag,np.log10(SNR),s=4,alpha=0.1,label='Stars')
plt.scatter(deep_mag,np.log10(deep_SNR),s=4,alpha=0.1,label='Deep')
plt.xlim(17,27)
plt.legend()
plt.savefig('zoomSNRcomp.png')
plt.clf()

plt.scatter(mag,SNR,s=4,alpha=0.1,label='Stars')
plt.scatter(deep_mag,deep_SNR,s=4,alpha=0.1,label='Deep')
plt.legend()
plt.savefig('newSNRcomp.png')
plt.clf()

SNe_median, SNe_median_edge, SNe_median_number = stats.binned_statistic(deep_mag,deep_SNR,statistic='median',bins=20,range=(20,30))
plt.hlines(bin_median[:-4], bin_median_edge[:-5], bin_median_edge[1:-4], colors='red',label='Stars')
plt.hlines(SNe_median, SNe_median_edge[:-1], SNe_median_edge[1:], colors='blue',label='SNANA')
#plt.plot(bin_median, bin_median_edge[:-1], label='Stars')
#plt.plot(SNe_median, SNe_median_edge[:-1], label='SNANA')
plt.legend()
plt.savefig('binnedSNRcomp.png')
plt.clf()


plt.hlines(SNe_median/bin_median, SNe_median_edge[:-1], SNe_median_edge[1:])
plt.savefig('ratioSNR.png')
plt.clf()
