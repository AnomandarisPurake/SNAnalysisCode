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
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import glob

no_of_scas = 18
snfile = fits.open('../truth/WFIRST_AKARI_SETEXP_KEVIN_HEAD.FITS')
sndata = snfile[1].data

def coadd_truth(filename):
    name = '../truth/paper_final_index_' + filename[-15:-4] + '.fits.gz'
    gals = fits.open(name)
    gdata = gals[1].data
    gdata = gdata[gdata['gal_star'] == 0]
    #print(len(data[data['mag'][:,3] != 0]))
    return gdata[gdata['mag'][:,3] != 0]

def truth_gals(ra_min, ra_max, dec_min, dec_max):
    data = sndata
    true_gals = data[(data['ra'] > ra_min) & (data['ra'] < ra_max) & (data['dec'] > dec_min) & (data['dec'] < dec_max)]
    return true_gals


galfile = fits.open('../truth/paper_final_truth_gal.fits')
galdata = galfile[1].data
sne = dict()
for f in glob.glob('paper_final_Y*.cat'):
    detected = 0
    cfile = f
    cdetections = pd.read_csv(cfile, skiprows=54, sep='\s+', header=None)
    cdetectx = np.array(cdetections[46])
    cdetecty = np.array(cdetections[47])
    cdetectra = np.array(cdetections[48])
    cdetectdec = np.array(cdetections[49])
    cdetectcoords = SkyCoord(cdetectra, cdetectdec, frame='fk5', unit=u.deg)
    ra_min = min(cdetectra) #* np.pi / 180
    ra_max = max(cdetectra) #* np.pi / 180
    dec_min = min(cdetectdec) #* np.pi / 180
    dec_max = max(cdetectdec) #* np.pi / 180
    data = truth_gals(ra_min, ra_max, dec_min, dec_max)
    detected = 0
    correct = 1
    

    old_correct = 1
        
    allgals = coadd_truth(f)
    allgalscoords = SkyCoord(allgals['ra'], allgals['dec'], frame='fk5', unit=u.deg)

    
    for i in range(len(data)):

        snid = data[i]['snid']
        if snid not in sne:
            gal = allgals[allgals['ind'] == data[i]['hostgal_objid']]
            if len(gal) != 0:
                gal_index = np.nonzero(allgals['ind'] == data[i]['hostgal_objid'])
                sncoords = SkyCoord(data[i]['ra'], data[i]['dec'], frame='fk5', unit=u.deg)
                galredshift = galdata[gal['ind'][0]]['z']
                all_indices, d2d, d3d = cdetectcoords.match_to_catalog_sky(allgalscoords)
                if gal_index in all_indices:
                    detected = 1
                

                #Old-style
                host_index, snd2d, snd3d = sncoords.match_to_catalog_sky(cdetectcoords)
                old_matchedgal = allgals[all_indices[host_index]]
                if old_matchedgal['ind'] == gal['ind']:
                    old_correct = 1
                    old_redshift = galredshift
                else:
                    old_correct = 0
                    old_redshift = galdata[old_matchedgal['ind']]['z']
                        
                sne[snid] = {'total': 1, 'correct': detected and correct, 'old': detected and old_correct, 'found': detected,
                            'redshift': data[i]['redshift_final'], 'galmag': gal['mag'][0], 'galredshift': galredshift, 'old_redshift': old_redshift}
        detected = 0
        correct = 1

        old_correct = 1

        #g = open(fname + str(pointing) + '_' + str(sca) + '.reg', mode='w')
        #for i in range(len(detectx)):
        #    g.write('circle ' + str(detectx[i]) + ' ' + str(detecty[i]) + ' 3.0 # color=red width = 3 \n')
        #g.close()
print(sne)
results = open('HostCoadd.txt', mode='w')
results.write(str(sne))
results.close()
found = np.array([])
notfound = np.array([])
correct = np.array([])
for v in sne.values():
    if float(v['old']) / float(v['total']) >= 0.5:
        correct = np.append(correct, v['galmag'])
    elif float(v['found']) / float(v['total']) >= 0.5:
        found = np.append(found, v['galmag'])
    else:
        notfound = np.append(notfound, v['galmag'])
bins = np.arange(19,27,0.5)
cbin = np.histogram(correct,bins)[0]
fbin = np.histogram(found,bins)[0]
nbin = np.histogram(notfound,bins)[0]
tbin = cbin + fbin + nbin
cfrac = cbin/tbin
ffrac = fbin/tbin
nfrac = nbin/tbin
plt.hist(bins[:-1],bins,weights=cfrac)
plt.ylim(0,1)
plt.ylabel('Correct fraction')
plt.xlabel('Y Band Magnitude')
plt.title('Correctly Identified Hosts')
plt.savefig('HostCoaddCorrect.png')
plt.clf()

cgrs = np.array([])
cwrs = np.array([])
fgrs = np.array([])
fwrs = np.array([])
nfgrs = np.array([])
nfwrs = np.array([])
for v in sne.values():
    if v['redshift'] < 10:
        if float(v['old']) / float(v['total']) >= 0.5:
            cgrs = np.append(cgrs, v['galredshift'])
            cwrs = np.append(cwrs, v['old_redshift'])
        elif float(v['found']) / float(v['total']) >= 0.5:
            fgrs = np.append(fgrs, v['galredshift'])
            fwrs = np.append(fwrs, v['old_redshift'])
        else:
            nfgrs = np.append(nfgrs, v['galredshift'])
            nfwrs = np.append(nfwrs, v['old_redshift'])

plt.scatter(cgrs,cwrs,label='Host correctly identified')
plt.scatter(fgrs,fwrs,label='Host detected but not correctly identified')
plt.scatter(nfgrs,nfwrs,label='Host not detected') 
plt.xlabel('True z')
plt.ylabel('z of Associated Galaxy')
plt.legend()
plt.tight_layout()
high = np.arange(0,3,0.1)
low = np.arange(0.2,3,0.1)
plt.plot(high, 0.2+high, color='red')
plt.plot(low, low - 0.2, color='red')
plt.gcf().subplots_adjust(bottom=0.15)



plt.savefig('OldHostCoaddRedshifts.png')


cgrs = np.array([])
nfgrs = np.array([])
fgrs = np.array([])
outlier = 0

for v in sne.values():
    if v['redshift'] < 10:
        if np.abs(v['galredshift'] - v['old_redshift']) > 0.2:
            outlier = outlier + 1
        elif float(v['old']) / float(v['total']) >= 0.5:
            cgrs = np.append(cgrs, v['galredshift'])
            cwrs = np.append(cwrs, v['old_redshift'])
        elif float(v['found']) / float(v['total']) >= 0.5:
            fgrs = np.append(fgrs, v['galredshift'])
            fwrs = np.append(fwrs, v['old_redshift'])
        else:
            nfgrs = np.append(nfgrs, v['galredshift'])
            nfwrs = np.append(nfwrs, v['old_redshift'])

print(len(cgrs))
print(len(fgrs))
print(len(nfgrs))
print(outlier)
