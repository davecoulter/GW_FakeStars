import glob
import numpy as np
from astropy.io import fits
import os
from txtobj import txtobj
import re
from runsex import runsex
import sys
import astropy.table as at

print("running SExtractor")


sextable = runsex(file_association.image_file, segmapname=file_association.segmap_file,
                                  zpt=fits.getval(file_association.image_file, 'ZPTMAG'))

image_id = file_association.image_file.split('/')[-1].split('.')[0]
glade_files = glob.glob('TileStats/*%s*txt'%image_id)
for gf in glade_files:
    #if len(glade_files):
    #	glade_file = glade_files[0]
    #else:
    #	raise RuntimeError('temporarily raising error here! can\'t find GLADE file!')
    good_ids,glade_ids,glade_bmags = np.array([]),np.array([]),np.array([])

    try:  glade = at.Table.read(gf,format='ascii.ecsv')
    except: continue

    for i in range(len(sextable.NUMBER)):
        sep = astCoords.calcAngSepDeg(glade['Galaxy_RA'],glade['Galaxy_Dec'],sextable.X_WORLD[i],sextable.Y_WORLD[i])*3600.
        if not len(sep):
            break
        if min(sep) < 2:
            good_ids = np.append(good_ids,sextable.NUMBER[i])
            glade_ids = np.append(glade_ids,glade['Galaxy_ID'][sep == np.min(sep)][0])
            glade_bmags = np.append(glade_bmags,glade['B'][sep == np.min(sep)][0])
if not len(glade_ids): continue