#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import wcs
import sys
import runsex


# In[ ]:


no_of_scas = 18
no_of_pointings = 44
indexname = '/hpc/group/cosmology/phy-lsst/kxw/akari_long/truth/akari_match_index_Y106_'
fname = '/hpc/group/cosmology/phy-lsst/kxw/akari_long/images/akari_match_Y106_'
sne = dict()
for sca in range(1, 1+no_of_scas):
    for pointing in range(no_of_pointings):
        file = fname + str(pointing) + '_' + str(sca) + '.fits'
        detections = runsex.runsex(file)
        detectx = detections.X_IMAGE
        detecty = detections.Y_IMAGE
        detected = False
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
                for j in range(len(detectx)):
                    xdisp = detectx[j] - snx
                    ydisp = detecty[j] - sny
                    if np.abs(xdisp) < np.abs(galxdisp)/2 and np.abs(ydisp) < np.abs(galydisp)/2:
                        detected = True
                sne[(snid, sca, pointing)] = {'detected': detected, 'mag': data[i]['mag']}
                detected = False
        a.close()
        b.close()
results = open('results.txt', mode='w')
results.write(str(sne))
results.close()
found = np.array([])
notfound = np.array([])
for v in sne.values():
    if v['detected']:
        found = np.append(found, v['mag'])
    else:
        notfound = np.append(notfound, v['mag'])
plt.hist(found, histtype = 'step')
plt.hist(notfound, histtype = 'step')
plt.savefig('Basic.png')
    

