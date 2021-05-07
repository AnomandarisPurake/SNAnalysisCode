from astropy.io import fits
import numpy as np
import sys
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord


# In[ ]:



indexname = '/hpc/group/cosmology/phy-lsst/kxw/akari_long/truth/akari_match_index_Y106_'
fname = '/hpc/group/cosmology/phy-lsst/kxw/akari_long/images/akari_match_Y106_'
cname = '/hpc/group/cosmology/phy-lsst/kxw/akari_long/images/akari_match_Y106_coadd_six_'

no_of_scas = 18
no_of_pointings = 44
for sca in range(1, 1+no_of_scas):
    for pointing in range(1,no_of_pointings):
        file = fname + str(pointing) + '_' + str(sca) + '.delme.txt'
        detections = pd.read_csv(file, skiprows=29, sep='\s+', header=None)
        detectx = np.array(detections[0])
        detecty = np.array(detections[1])
        
        g = open(fname + str(pointing) + '_' + str(sca) + '.reg', mode='w')
        for i in range(len(detectx)):
            g.write('circle ' + str(detectx[i]) + ' ' + str(detecty[i]) + ' 4.0 # color=red width = 2 \n')
        g.close()

no_of_scas = 18
no_of_pointings = 1
for sca in range(2, 1+no_of_scas):
    for pointing in range(1,no_of_pointings+1):
        file = fname + str(pointing) + '_' + str(sca) + '.delme.txt'
        detections = pd.read_csv(file, skiprows=29, sep='\s+', header=None)
        detectx = np.array(detections[12])
        detecty = np.array(detections[13])
        detectcoords = SkyCoord(detectx, detecty, frame='fk5', unit=u.deg)
        
        cfile = cname + str(sca) + '_drz.delme.txt'
        cdetections = pd.read_csv(cfile, skiprows=29, sep='\s+', header=None)
        cdetectx = np.array(cdetections[0])
        cdetecty = np.array(cdetections[1])
        cdetectra = np.array(cdetections[12])
        cdetectdec = np.array(cdetections[13])
        cdetectcoords = SkyCoord(cdetectra, cdetectdec, frame='fk5', unit=u.deg)
        g = open(cname + str(sca) + '.reg', mode='w')
        for i in range(len(cdetectx)):
            seps = cdetectcoords[i].separation(detectcoords).arcsec
            sepssort = np.argsort(seps)
            if seps[sepssort[0]] < 0.16:
                g.write('circle ' + str(cdetectx[i]) + ' ' + str(cdetecty[i]) + ' 8.0 # color=red width = 2 \n')
            else:
                g.write('circle ' + str(cdetectx[i]) + ' ' + str(cdetecty[i]) + ' 8.0 # color=cyan width = 2 \n')
        g.close()
