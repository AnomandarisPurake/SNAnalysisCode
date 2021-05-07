from astropy.io import fits
import numpy as np
import sys
import runsex


# In[ ]:


no_of_scas = 18
no_of_pointings = 1
indexname = '/hpc/group/cosmology/phy-lsst/kxw/akari_long/truth/akari_match_index_Y106_'
fname = '/hpc/group/cosmology/phy-lsst/kxw/akari_long/images/akari_match_Y106_coadd_six_'
gals = dict()
for sca in range(1, 1+no_of_scas):
    for pointing in range(no_of_pointings):
        file = fname + str(sca) + '_drz.fits'
        detections = runsex.runsex(file)
        detectx = np.array(detections.X_IMAGE)
        detecty = np.array(detections.Y_IMAGE)
        g = open(fname + str(sca) + '.reg', mode='w')
        for i in range(len(detectx)):
            g.write('circle ' + str(detectx[i]) + ' ' + str(detecty[i]) + ' 6.0 # color=red width = 3 \n')
        g.close()
