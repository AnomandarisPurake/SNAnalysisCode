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
snfile = fits.open('/hpc/group/cosmology/phy-lsst/kxw/Matched_Gal_HEAD.FITS')
sndata = snfile[1].data
sne = dict()
for sca in range(1, 1+no_of_scas):
    for pointing in range(no_of_pointings):
        file = fname + str(pointing) + '_' + str(sca) + '.fits'
        detections = runsex.runsex(file)
        detectx = detections.X_IMAGE
        detecty = detections.Y_IMAGE
        detected = 0
        correct = 1
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
                    if (xdisp - galxdisp)**2 + (ydisp - galydisp)**2 < 3:
                        detected = 1
                    if xdisp**2 + ydisp**2 < galxdisp**2 + galydisp**2 - 3 and xdisp**2 + ydisp**2 > 3:
                        correct = 0
                if snid not in sne:
                    thissn = sndata[np.asarray(sndata['snid']).astype(int) == data[i]['ind']]
                    sne[snid] = {'total': 1, 'correct': detected and correct, 'found': detected, 'redshift': thissn['redshift_final'][0], 'galmag': gal['mag'][0]}
                else:
                    sne[snid]['total'] = sne[snid]['total'] + 1
                    sne[snid]['correct'] = sne[snid]['correct'] + (detected and correct)
                    sne[snid]['found'] = sne[snid]['found'] + detected
                detected = 0
                correct = 1
        a.close()
        b.close()
        #g = open(fname + str(pointing) + '_' + str(sca) + '.reg', mode='w')
        #for i in range(len(detectx)):
        #    g.write('circle ' + str(detectx[i]) + ' ' + str(detecty[i]) + ' 3.0 # color=red width = 3 \n')
        #g.close()
results = open('HostResults.txt', mode='a')
results.write(str(sne))
results.close()
found = np.array([])
notfound = np.array([])
correct = np.array([])
notcorrect = np.array([])
for v in sne.values():
    if float(v['found']) / float(v['total']) >= 0.5:
        found = np.append(found, v['redshift'])
    else:
        notfound = np.append(notfound, v['redshift'])
plt.hist(found, histtype = 'step')
plt.hist(notfound, histtype = 'step')
plt.savefig('HostFound.png')
plt.clf()
for v in sne.values():
    if float(v['correct']) / float(v['total']) >= 0.5:
        correct = np.append(correct, v['redshift'])
    else:
        notcorrect = np.append(notcorrect, v['redshift'])
plt.hist(correct, histtype = 'step')
plt.hist(notcorrect, histtype = 'step')
plt.savefig('HostCorrect.png')

    

