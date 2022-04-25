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
import plotsetup
import matplotlib.gridspec as gridspec

# In[ ]:


no_of_scas = 18
no_of_pointings = 44
indexname = '../truth/paper_final_index_Y106_'
fname = 'paper_final_Y106_'
stars = dict()
found = np.array([])
notfound = np.array([])


stars = eval(open('StarSNRResults.txt',mode='r').read())

SNR = np.array([])
mag = np.array([])
Smag = np.array([])
for v in stars.values():
    if v['found'] >= 10:
        SNR = np.append(SNR, v['SNR']/v['found'])
        mag = np.append(mag, v['mag'])
        Smag = np.append(Smag, v['Smag'])

plotsetup.halfpaperfig()

bin_median, bin_median_edge, bin_median_number = stats.binned_statistic(mag,SNR,statistic='median',bins=20,range=(20,30))
bin_std, bin_std_edge, bin_std_number = stats.binned_statistic(mag,SNR,statistic='std',bins=20,range=(20,30))
    
other = pd.read_csv('OUT_FOR_KEVIN.OUTLIER.TEXT',skiprows=15, sep='\s+',header=None)
logSNR = np.array(other[17])
fluxcal = np.array(other[15])
band = np.array(other[5])
field = np.array(other[3])
deep = np.intersect1d(np.asarray(band == 'Y').nonzero()[0],np.asarray(field == 'DEEP').nonzero()[0])
deep_mag = 27.5 - 2.5 * np.log10(fluxcal[deep])
deep_SNR = np.power(10,logSNR[deep])
plt.scatter(mag,SNR,s=4,alpha=0.2, color='blue')
plt.scatter(deep_mag,deep_SNR,s=4,alpha=0.1, color='orange')
plt.scatter([],[], label='Stars', color='blue')
plt.scatter([],[], label='SNANA', color='orange')

SNe_median, SNe_median_edge, SNe_median_number = stats.binned_statistic(deep_mag,deep_SNR,statistic='median',bins=20,range=(20,30))
SNe_std, SNe_std_edge, SNe_std_number = stats.binned_statistic(deep_mag,deep_SNR,statistic='std',bins=20,range=(20,30))
#plt.hlines(bin_median[:-4], bin_median_edge[:-5], bin_median_edge[1:-4], colors='red',label='Stars')
#plt.hlines(SNe_median, SNe_median_edge[:-1], SNe_median_edge[1:], colors='blue',label='SNANA')
#plt.plot(bin_median, bin_median_edge[:-1], label='Stars')
#plt.plot(SNe_median, SNe_median_edge[:-1], label='SNANA')
star_y = np.ravel(list(zip(bin_median,bin_median)))
star_x = np.ravel(list(zip(bin_median_edge[:-1],bin_median_edge[:-1] + 0.5)))
SN_y = np.ravel(list(zip(SNe_median,SNe_median)))
SN_x = np.ravel(list(zip(SNe_median_edge[:-1],SNe_median_edge[:-1] + 0.5)))
plt.plot(star_x, star_y, label='Stars Median', color='black')
plt.plot(SN_x, SN_y, label='SNANA Median', color='red')
#plt.errorbar(bin_mean_edge[:-1] + 0.5, bin_mean, yerr = bin_std, fmt='none', label='Star error')
#plt.errorbar(bin_mean_edge[:-1] + 0.5, SNe_mean, yerr = SNe_std, fmt='none', label='SNANA error')
plt.legend()
plt.xlabel('Y Band Magnitude')
plt.ylabel('S/N')
plt.xlim(20,26.5)
plt.ylim(1,300)
plt.yscale('log')
plt.tight_layout()
plt.gcf().subplots_adjust(bottom=0.15)
plt.savefig('binnedSNRcomp.png')
plt.clf()
