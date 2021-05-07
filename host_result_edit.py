#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import wcs
import sys


# In[ ]:


no_of_scas = 18
no_of_pointings = 44
indexname = '/hpc/group/cosmology/phy-lsst/kxw/akari_long/truth/akari_match_index_Y106_'
fname = '/hpc/group/cosmology/phy-lsst/kxw/akari_long/images/akari_match_Y106_'
snfile = fits.open('/hpc/group/cosmology/phy-lsst/kxw/Matched_Gal_HEAD.FITS')
sndata = snfile[1].data
a = open('HostResults.txt',mode='r').read()
sne = eval(a)
for sca in range(1, 1+no_of_scas):
    for pointing in range(no_of_pointings):
        index = indexname + str(pointing) + '_' + str(sca) + '_sn.fits'
        galindex = indexname + '0_'+ str(sca) + '.fits'
        a = fits.open(index)
        b = fits.open(galindex)
        data = a[1].data
        for i in range(len(data)):
            snid = data[i]['ind']
            if data[i]['mag'] < 30:
                gal = b[1].data[b[1].data['ind'] == data[i]['hostid']]
                snx = data[i]['x']
                sny = data[i]['y']
                galxdisp = gal['x'] - snx
                galydisp = gal['y'] - sny
                if len(gal) != 0:
                    sne[snid]['galmag'] = gal['mag'][0]
                else:
                    sne[snid]['galmag'] = sndata[np.asarray(sndata['snid']).astype(int) == snid]['hostgal_mag_y'][0]
        a.close()
        b.close()
        #g = open(fname + str(pointing) + '_' + str(sca) + '.reg', mode='w')
        #for i in range(len(detectx)):
        #    g.write('circle ' + str(detectx[i]) + ' ' + str(detecty[i]) + ' 3.0 # color=red width = 3 \n')
        #g.close()
results = open('HostResults.txt', mode='w')
results.write(str(sne))
results.close()
found = np.array([])
notfound = np.array([])
correct = np.array([])
for v in sne.values():
    if float(v['correct']) / float(v['total']) >= 0.5:
        correct = np.append(correct, v['galmag'])
    elif float(v['found']) / float(v['total']) >= 0.5:
        found = np.append(found, v['galmag'])
    else:
        notfound = np.append(notfound, v['galmag'])
bins = np.arange(17,26,1)
cbin = np.histogram(correct,np.arange(17,26,1))[0]
fbin = np.histogram(found,np.arange(17,26,1))[0]
nbin = np.histogram(notfound,np.arange(17,26,1))[0]
tbin = cbin + fbin + nbin
np.seterr(divide='ignore')
cfrac = cbin/tbin
ffrac = fbin/tbin
nfrac = nbin/tbin
print(cfrac)
plt.hist(bins[:-1],bins,weights=cfrac)
plt.ylim(0,1)
plt.title('Correctly Identified Hosts')
plt.savefig('Correct.png')
plt.clf()
plt.hist(bins[:-1],bins,weights=ffrac)
plt.ylim(0,1)
plt.title('Found but Not Closest Hosts')
plt.savefig('Found.png')
plt.clf()
plt.hist(bins[:-1],bins,weights=nfrac)
plt.ylim(0,1)
plt.title('Not Found Hosts')
plt.savefig('NotFound.png')

    

