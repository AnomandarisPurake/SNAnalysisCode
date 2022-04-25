#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from astropy.io import fits
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import pylab
import sys
import plotsetup

# In[ ]:

plotsetup.halfpaperfig()
a = fits.open('Paper_KXW_HEAD.FITS')
b = fits.open('paper_final_truth_gal.fits')
hdata = a[1].data
adata = b[1].data
hostmags = np.array([])
for i in range(len(hdata)):
    host = adata[adata['gind'] == hdata['hostgal_objid'][i]]
    if host['z'] < 1:
        hostmags = np.append(hostmags,host['Y106'])
plotsetup.halfpaperfig()
plt.hist(adata['Y106'],bins=30)
plt.title('Y Band Mags of All Galaxies')
plt.xlabel('Magnitude')
plt.ylabel('Number of Galaxies')
plt.tight_layout()
plt.gcf().subplots_adjust(bottom=0.15)
plt.savefig('GalMagsY.png')
plt.clf()

plotsetup.halfpaperfig()
plt.hist(hostmags,bins=30)
plt.title('Y Band Mags of Host Galaxies')
plt.xlabel('Magnitude')
plt.ylabel('Number of Galaxies')
plt.tight_layout()
plt.gcf().subplots_adjust(bottom=0.15)
plt.savefig('HostMagsY.png')
plt.clf()

plotsetup.halfpaperfig()
plt.hist(adata['Y106'],bins=30,range=(15.4,29.5),density=True,histtype='step',label='All Galaxies')
plt.hist(hostmags,bins=30,range=(15.4,29.5),density=True,histtype='step',label='Host Galaxies')
plt.title('Galaxies Y Band Magnitudes Comparison')
plt.xlabel('Magnitude')
plt.ylabel('Fraction of Galaxies')
plt.legend()
plt.tight_layout()
plt.gcf().subplots_adjust(bottom=0.15)
plt.savefig('MagsY.png')
plt.clf()

plotsetup.halfpaperfig()
plt.hist(adata['z'],bins=30, range=(0,5), density=True, histtype='step', label='Galaxies')
plt.hist(a[1].data['redshift_final'],bins=30, range=(0,5), density=True, histtype='step', label='Supernovae')
plt.title('Redshift distribution of galaxies')
plt.legend()
plt.xlabel('z')
plt.ylabel('Number of objects')
plt.tight_layout()
plt.gcf().subplots_adjust(bottom=0.15)
plt.savefig('GalSNZ.png')
