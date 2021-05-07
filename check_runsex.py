import runsex
import pandas as pd
import numpy as np

cnames = ('X_IMAGE','Y_IMAGE','MAG_AUTO','MAGERR_AUTO','FLAGS','BACKGROUND','FWHM_IMAGE','THETA_IMAGE','ELONGATION','CLASS_STAR','FLUX_AUTO','FLUXERR_AUTO','X_WORLD','Y_WORLD','A_IMAGE','B_IMAGE','A_WORLD','B_WORLD','FWHM_WORLD','MU_MAX','FLUX_MAX')
columns = np.arange(21)

text = pd.read_csv('fiducial_Y106_2_1.delme.txt',sep='\s+',names=cnames,usecols=columns)
a = runsex.runsex('fiducial_Y106_2_1.fits')
print(np.array_equal(a.X_IMAGE,np.array(text.X_IMAGE)))
print(np.array_equal(a.X_WORLD,np.array(text.X_WORLD)))
print(np.array_equal(a.A_IMAGE,np.array(text.A_IMAGE)))
print(a.A_IMAGE)
print(text.A_IMAGE)
