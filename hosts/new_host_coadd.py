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

def Cij(sextractor_results):
    # A_WORLD,B_WORLD,THETA_IMAGE
    a,b,theta = np.array(sextractor_results[54]),np.array(sextractor_results[55]),np.array(sextractor_results[53]) 
    theta *= np.pi/180
    # Cij ~ /pixel^2
    Cxx = np.cos(theta)**2/a**2 + np.sin(theta)**2/b**2
    Cyy = np.sin(theta)**2/a**2 + np.cos(theta)**2/b**2
    Cxy = 2*np.cos(theta)*np.sin(theta)*(1/a**2+1/b**2)
        
    try:
        assert(sextractor_results['Cxx'] is not None)
    except:
        print("Adding Cij columns to sextractor results")
        sextractor_results['Cxx'] = Cxx
        sextractor_results['Cyy'] = Cyy
        sextractor_results['Cxy'] = Cxy
    
    return Cxx,Cyy,Cxy

def sep(sn,sextractor_results):
    # pix
    xg = np.array(sextractor_results[48])
    yg = np.array(sextractor_results[49])
    allgals = SkyCoord(xg, yg, frame='fk5', unit=u.deg)
    seps = []
    for i in range(len(sn)):
        # pix
        sncoords = SkyCoord(sn[i]['ra'], sn[i]['dec'], frame='fk5', unit=u.deg)
        dists = sncoords.separation(allgals).arcsec
        angles = sncoords.position_angle(allgals).to(u.rad)
        xoffsets = dists * np.sin(angles)
        yoffsets = dists * np.cos(angles)
        seps.append([xoffsets,yoffsets])
    return seps

def R(sextractor_results,sep):
    """
    Determine a host to SN using galaxies elliptical shapes determined using SExtractor 
    # Sullivan https://arxiv.org/pdf/astro-ph/0605455.pdf
    # R^2 = Cxx xr^2 + Cyy yr^2 + Cxy xr*yr [dimensionless, units of galaxy ellipse size]
    # Cij ~ /[pixel^2] ~ f(a,b,theta)
    # xr/yr ~ [pixel] ~ pix_sn - pix_gal
    """

    xr,yr = sep
    Cxx = sextractor_results['Cxx']
    Cyy = sextractor_results['Cyy']
    Cxy = sextractor_results['Cxy']
    
    Reff = Cxx*xr**2 + Cyy*yr**2 + Cxy*xr*yr
    
    Reff.name = "R"
    Reff.unit = "" # dimensionless ~ units of the galaxy ellipse
    
##    try:
##        assert(sextractor_results['R'] is not None)
##        # After 0th SN need to replace the R column to match up gals to this ith SN
##        print("Replacing R column in sextractor_results")
##        sextractor_results.replace_column("R",Reff)
##    except:
##        # the 0th SN wont have the R column
##        print("Adding R column to sextractor results")
##        sextractor_results.add_column(Reff)
    
    host_idx = np.argmin(Reff)
##    host = sextractor_results[host_idx]
    
    return host_idx

##if __name__ == "__main__":
##
##    # dither/sca ~ epoch/detector ~ 44*18 ~ 792
##    dithers = [i for i in range(44)] # [0,43]
##    scas = [i for i in range(1,19)] # [1,18]
##
##    for dither in dithers:
##        for sca in scas:
##            
##            # get the data/parse into needed objects
##            print("dither ~ {}, sca ~ {}".format(dither,sca))
##            sex_results, ass_truths = data_dither_sca(dither,sca)
##            sex_ex,coadd_ex = sex_results
##            ass_ex,sn_ex,star_ex = ass_truths
##            sn_data = sn_ex[1].data
##
##            # each dither_sca has multiple planted SN
##            print("{} planted SNe".format(len(sn_data)))
##            
##            # find hosts using coadded images or the single dither_sca image
##            use_coadds = True
##            if use_coadds:
##                sextractor = coadd_ex
##            else:
##                sextractor = sex_ex
##
##            # determine SN host galaxies, minimum R (Sullivan https://arxiv.org/pdf/astro-ph/0605455.pdf)
##            Cij(sextractor)
##            seps = sep(sn_data,sextractor) # listlike for N SNe each a table of the diff from sextractor gals 
##            print("{} seps".format(len(seps)))
##            for i in range(len(sn_data)):
##                print("sn {}".format(i))
##                sepi = seps[i]
##                hosti = R(sextractor,sepi)  
##                
##                print("host",i,"R ~ ",hosti["R"])
##                if use_coadds:
##                    saveas = f'{dither}_{sca}_coadded_host{i}.txt'
##                else:
##                    saveas = f'{dither}_{sca}_host{i}.txt'
##
##                ascii.write(hosti, saveas, overwrite=True)

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
        
    allgals = coadd_truth(f)
    allgalscoords = SkyCoord(allgals['ra'], allgals['dec'], frame='fk5', unit=u.deg)

    # Steve stuff
    # determine SN host galaxies, minimum R (Sullivan https://arxiv.org/pdf/astro-ph/0605455.pdf)
    Cij(cdetections)
    seps = sep(data,cdetections) # listlike for N SNe each a table of the diff from sextractor gals

    for i in range(len(data)):
        sepi = seps[i]
        hosti = R(cdetections,sepi)
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
                    
                matchedgal = allgals[all_indices[hosti]]
                if matchedgal['ind'] == gal['ind']:
                    correct = 1
                    wrongredshift = galredshift
                        
                else:
                    correct = 0
                    wrongredshift = galdata[matchedgal['ind']]['z']
                        
                sne[snid] = {'total': 1, 'correct': detected and correct, 'found': detected,
                            'redshift': data[i]['redshift_final'], 'galmag': gal['mag'][0], 'galredshift': galredshift, 'wrongredshift': wrongredshift}
        detected = 0
        correct = 1
            
        #g = open(fname + str(pointing) + '_' + str(sca) + '.reg', mode='w')
        #for i in range(len(detectx)):
        #    g.write('circle ' + str(detectx[i]) + ' ' + str(detecty[i]) + ' 3.0 # color=red width = 3 \n')
        #g.close()
results = open('NewHostCoadd.txt', mode='a')
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
bins = np.arange(17,27,1)
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
plt.savefig('NewHostCoaddCorrect.png')
plt.clf()
plt.hist(bins[:-1],bins,weights=ffrac)
plt.ylim(0,1)
plt.title('Found but Not Closest Hosts')
plt.savefig('NewHostCoaddFound.png')
plt.clf()
plt.hist(bins[:-1],bins,weights=nfrac)
plt.ylim(0,1)
plt.title('Not Found Hosts')
plt.savefig('NewHostCoaddNotFound.png')
plt.clf()
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
plt.ylabel('z of Associated Galaxy')
plt.legend()
plt.tight_layout()
plt.gcf().subplots_adjust(bottom=0.15)
plt.savefig('NewHostCoaddRedshifts.png')
plt.clf()

print(len(cgrs))
print(len(fgrs))
print(len(nfgrs))    

