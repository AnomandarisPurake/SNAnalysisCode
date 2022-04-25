#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import sys
import plotsetup
import matplotlib.gridspec as gridspec

# In[ ]:
allgals = fits.open('paper_final_truth_gal.fits')

#Single Image General Galaxy Results
plotsetup.halfpaperfig()
fig, host = plt.subplots()
#second = host.twinx()

gals = eval(open('SingleStarResults.txt',mode='r').read())
found = np.array([])
notfound = np.array([])
for v in gals.values():
    if float(v['found']) / float(v['total']) >= 0.5:
        found = np.append(found, v['galmag'])
    else:
        notfound = np.append(notfound, v['galmag'])
bins = np.arange(22,28,0.1)
fbin = np.histogram(found,bins)[0]
nbin = np.histogram(notfound,bins)[0]
tbin = fbin + nbin
ffrac = fbin/tbin
nfrac = nbin/tbin
host.hist(bins[:-1],bins,weights=ffrac,histtype='step',label='Single Images')
host.set_ylim(0,1)

#Coadd General Galaxy Results
gals = eval(open('CoaddStarResults.txt',mode='r').read())
found = np.array([])
notfound = np.array([])
for v in gals.values():
    if float(v['found']) / float(v['total']) >= 0.5:
        found = np.append(found, v['galmag'])
    else:
        notfound = np.append(notfound, v['galmag'])
bins = np.arange(22,28,0.1)
fbin = np.histogram(found,bins)[0]
nbin = np.histogram(notfound,bins)[0]
tbin = fbin + nbin
ffrac = fbin/tbin
nfrac = nbin/tbin
host.hist(bins[:-1],bins,weights=ffrac,histtype='step',label='Coadded Images')
host.legend(loc=2)
#second.hist(allgals[1].data['Y106'],bins,alpha=0.4,color='0.8') 
host.set_xlabel('Magnitude')
host.set_ylabel('Detection efficiency')
#second.set_ylabel('Number of Stars')
plt.minorticks_on()
plt.tight_layout()
plt.gcf().subplots_adjust(bottom=0.15)
plt.savefig('StarEfficiency.png')
plt.clf()
