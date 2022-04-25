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
second = host.twinx()

gals = eval(open('SingleResults.txt',mode='r').read())
found = np.array([])
notfound = np.array([])
for v in gals.values():
    if float(v['found']) / float(v['total']) >= 0.5:
        found = np.append(found, v['galmag'])
    else:
        notfound = np.append(notfound, v['galmag'])
bins = np.arange(22,30,0.1)
fbin = np.histogram(found,bins)[0]
nbin = np.histogram(notfound,bins)[0]
tbin = fbin + nbin
ffrac = fbin/tbin
nfrac = nbin/tbin
host.hist(bins[:-1],bins,weights=ffrac,histtype='step',label='Single Images')
host.set_ylim(0,1)

#Coadd General Galaxy Results
gals = eval(open('CoaddResults.txt',mode='r').read())
found = np.array([])
notfound = np.array([])
for v in gals.values():
    if float(v['found']) / float(v['total']) >= 0.5:
        found = np.append(found, v['galmag'])
    else:
        notfound = np.append(notfound, v['galmag'])
bins = np.arange(22,30,0.1)
fbin = np.histogram(found,bins)[0]
nbin = np.histogram(notfound,bins)[0]
tbin = fbin + nbin
ffrac = fbin/tbin
nfrac = nbin/tbin
host.hist(bins[:-1],bins,weights=ffrac,histtype='step',label='Coadded Images')
host.legend()
second.hist(allgals[1].data['Y106'],bins,alpha=0.4,color='0.8') 
host.set_xlabel('Magnitude')
host.set_ylabel('Detection efficiency')
second.set_ylabel('Number of galaxies')
plt.minorticks_on()
plt.tight_layout()
plt.gcf().subplots_adjust(bottom=0.15)
plt.savefig('GalaxyEfficiency.png')
plt.clf()


sne = eval(open('HostCoadd.txt',mode='r').read())    
cgrs = np.array([])
cwrs = np.array([])
fgrs = np.array([])
fwrs = np.array([])
nfgrs = np.array([])
nfwrs = np.array([])

for v in sne.values():
    if float(v['correct']) / float(v['total']) >= 0.5:
        cgrs = np.append(cgrs, v['galredshift'])
        cwrs = np.append(cwrs, v['wrongredshift'])
    elif float(v['found']) / float(v['total']) >= 0.5:
        fgrs = np.append(fgrs, v['galredshift'])
        fwrs = np.append(fwrs, v['wrongredshift'])
    else:
        nfgrs = np.append(nfgrs, v['galredshift'])
        nfwrs = np.append(nfwrs, v['wrongredshift'])

plt.scatter(cgrs,cwrs,label='Host correctly identified')
plt.scatter(fgrs,fwrs,label='Host detected but not correctly identified')
plt.scatter(nfgrs,nfwrs,label='Host not detected') 
plt.xlabel('True z')
plt.ylabel('Calculated z')
plt.legend()
plt.tight_layout()
plt.gcf().subplots_adjust(bottom=0.15)
plt.savefig('HostCoaddRedshifts.png')
