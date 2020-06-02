import glob
import numpy as np
from astropy.io import fits
import os
from txtobj import txtobj
import re
from runsex import runsex
import sys
import astropy.table as at
from astLib import astCoords

print("running SExtractor")


image_file = "s005aae0078.i.ut190425.1083_stch_1.sw.fits"
sextable = runsex(image_file,
                  segmapname=image_file.replace('.fits', '.check.fits'),
                  zpt=fits.getval(image_file, 'ZPTMAG'))

print(sextable)

field_name = image_file.split('.')[0]
glade_files = glob.glob('*%s*txt' % field_name)

for gf in glade_files:
    good_ids, glade_ids, glade_bmags = np.array([]),np.array([]),np.array([])

    try:
        glade = at.Table.read(gf, format='ascii.ecsv')
    except:
        continue

    for i in range(len(sextable.NUMBER)):

        # separation in arc sec
        sep = astCoords.calcAngSepDeg(glade['Galaxy_RA'],
                                      glade['Galaxy_Dec'],
                                      sextable.X_WORLD[i],
                                      sextable.Y_WORLD[i])*3600.

        # check if it's found in the sextable...
        if not len(sep):
            break

        # if within 2 arc sec
        if min(sep) < 2:
            good_ids = np.append(good_ids, sextable.NUMBER[i])
            glade_ids = np.append(glade_ids, glade['Galaxy_ID'][sep == np.min(sep)][0])
            glade_bmags = np.append(glade_bmags, glade['B'][sep == np.min(sep)][0])

    print(good_ids)
    print(glade_ids)
    print(glade_bmags)