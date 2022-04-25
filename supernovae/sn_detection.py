#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from astropy.io import fits
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from astropy.wcs import wcs
import sys
import runsex
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import plotsetup
import matplotlib.gridspec as gridspec

# In[ ]:

plotsetup.halfpaperfig()

found = eval(open('SNFound.txt',mode='r').read())
notfound = eval(open('SNNotFound.txt',mode='r').read())
plt.hist((found,np.append(found,notfound)), histtype = 'step',stacked=False,label=('SN detections','SN injections'),bins=30)
plt.xlabel('Supernova magnitude (Y band)')
plt.ylabel('Number of SNe')
plt.legend()
plt.tight_layout()
plt.gcf().subplots_adjust(bottom=0.15)
plt.savefig('SNFound.png')
plt.clf()
