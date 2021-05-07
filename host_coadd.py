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
import runsex
import astropy.units as u
from astropy.coordinates import SkyCoord

# In[ ]:


no_of_scas = 18
no_of_pointings = 44
indexname = '/hpc/group/cosmology/phy-lsst/kxw/akari_long/truth/akari_match_index_Y106_'
fname = '/hpc/group/cosmology/phy-lsst/kxw/akari_long/images/akari_match_Y106_coadd_six_'
snfile = fits.open('/hpc/group/cosmology/phy-lsst/kxw/Matched_Gal_HEAD.FITS')
sndata = snfile[1].data
galfile = fits.open('/hpc/group/cosmology/phy-lsst/kxw/akari_long/truth/akari_match_truth_gal.fits')
galdata = galfile[1].data
sne = dict()
badgals = 0
for sca in range(1, 1+no_of_scas):
    file = fname + str(sca) + '_drz.fits'
    detections = runsex.runsex(file)
    detectx = np.array(detections.X_WORLD)
    detecty = np.array(detections.Y_WORLD)
    detectcoords = SkyCoord(detectx, detecty, frame='fk5', unit=u.deg)
    for pointing in range(no_of_pointings):
        detected = 0
        correct = 1
        index = indexname + str(pointing) + '_' + str(sca) + '_sn.fits'
        galindex = indexname + '0_' + str(sca) + '.fits'
        a = fits.open(index)
        b = fits.open(galindex)
        data = a[1].data
        allgals = b[1].data
        allgalscoords = SkyCoord(allgals['ra'], allgals['dec'], frame='fk5', unit=u.rad)
        for i in range(len(data)):
            snid = data[i]['ind']
            if snid not in sne:
                gal = b[1].data[b[1].data['ind'] == data[i]['hostid']]
                if len(gal) == 0:
                    badgals += 1
                    print(badgals)
                if len(gal) != 0:
                    sncoords = SkyCoord(data[i]['ra'], data[i]['dec'], frame='fk5', unit=u.rad)
                    galcoords = SkyCoord(gal['ra'], gal['dec'], frame='fk5', unit=u.rad)
                    dists = galcoords.separation(detectcoords).arcsec
                    distsort = np.argsort(dists)
                    if dists[distsort[0]] < 0.3:
                        detected = 1
                    sndists = sncoords.separation(detectcoords).arcsec
                    snsort = np.argsort(sndists)
                    if sndists[snsort[0]] < 0.1:
                        if not(snsort[0] == distsort[0] or snsort[1] == distsort[0]) or detected == 0:
                            matcher = detectcoords[snsort[1]]
                            alldists = matcher.separation(allgalscoords).arcsec
                            matchedgal = allgals[np.argsort(alldists)[0]]
                            if matchedgal['ind'] != gal['ind'][0]:
                                correct = 0
                    else:
                        if not(snsort[0] == distsort[0]) or detected == 0:
                            matcher = detectcoords[snsort[0]]
                            alldists = matcher.separation(allgalscoords).arcsec
                            matchedgal = allgals[np.argsort(alldists)[0]]
                            if matchedgal['ind'] != gal['ind'][0]:
                                correct = 0
                    galredshift = galdata[gal['ind'][0]]['z']
                    if correct == 0:
                        wrongredshift = galdata[matchedgal['ind']]['z']
                    else:
                        wrongredshift = galredshift
                    thissn = sndata[np.asarray(sndata['snid']).astype(int) == data[i]['ind']]
                    sne[snid] = {'total': 1, 'correct': detected and correct, 'found': detected, 'redshift': thissn['redshift_final'][0], 'galmag': gal['mag'][0], 'galredshift': galredshift, 'wrongredshift': wrongredshift}
            detected = 0
            correct = 1
        a.close()
        b.close()
        #g = open(fname + str(pointing) + '_' + str(sca) + '.reg', mode='w')
        #for i in range(len(detectx)):
        #    g.write('circle ' + str(detectx[i]) + ' ' + str(detecty[i]) + ' 3.0 # color=red width = 3 \n')
        #g.close()
results = open('HostCoadd.txt', mode='w')
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
cbin = np.histogram(correct,bins)[0]
fbin = np.histogram(found,bins)[0]
nbin = np.histogram(notfound,bins)[0]
tbin = cbin + fbin + nbin
cfrac = cbin/tbin
ffrac = fbin/tbin
nfrac = nbin/tbin
plt.hist(bins[:-1],bins,weights=cfrac)
plt.ylim(0,1)
plt.title('Correctly Identified Hosts')
plt.savefig('HostCoaddCorrect.png')
plt.clf()
plt.hist(bins[:-1],bins,weights=ffrac)
plt.ylim(0,1)
plt.title('Found but Not Closest Hosts')
plt.savefig('HostCoaddFound.png')
plt.clf()
plt.hist(bins[:-1],bins,weights=nfrac)
plt.ylim(0,1)
plt.title('Not Found Hosts')
plt.savefig('HostCoaddNotFound.png')
plt.clf()
grs = np.array([])
wrs = np.array([])
for v in sne.values():
    grs = np.append(grs, v['galredshift'])
    wrs = np.append(wrs, v['wrongredshift'])
plt.scatter(grs,wrs)
plt.xlabel('Redshift of host galaxy')
plt.ylabel('Redshift of closest galaxy')
plt.savefig('HostCoaddRedshifts.png')

    

