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
import astropy.coordinates as coord
from astropy import units as u
import time
from astropy.table import Table

def write_good_sexcat_ids(glade_file, image_file, good_ids, glade_ids, glade_bmags):

    ascii_ecsv_fname = "%s_sexcat_good.txt" % glade_file
    ascii_ecsv_fpath = "%s/%s" % (".", ascii_ecsv_fname)
    print("Creating `%s`" % ascii_ecsv_fpath)

    # Build ascii.ecsv formatted output
    cols = ['image_file', 'sexcat_id', 'glade_id', 'glade_B']
    dtype = ['U64', 'i4', 'i4', 'f8']
    result_table = Table(dtype=dtype, names=cols)
    meta = ["{key}={value}".format(key="glade_file", value=glade_file)]
    result_table.meta['comment'] = meta

    for sexcat_id, glade_id, b in zip(good_ids, glade_ids, glade_bmags):
        result_table.add_row([image_file, sexcat_id, glade_id, b])
    result_table.write(ascii_ecsv_fpath, overwrite=True, format='ascii.ecsv')

t1 = time.time()

print("running SExtractor")


dcmp_file = "s005aae0078.i.ut190425.1083_stch_1.sw.dcmp"
dcmp_header = fits.getheader(dcmp_file)
zpt = dcmp_header['ZPTMAG']

mask_file = " s005aae0078.i.ut190425.1083_stch_1.sw.mask.fits.gz"
mask_hdu = fits.open(mask_file)
mask_data = mask_hdu[0].data.astype('float')

image_file = "s005aae0078.i.ut190425.1083_stch_1.sw.fits"
field_name = image_file.split('.')[0]
glade_files = glob.glob('*%s*txt' % field_name)

# RUN SEXTRACTOR - this outputs the file *.sexcat
segmap_file = image_file.replace('.fits', '.check.fits')
sextable = runsex(image_file, segmapname=segmap_file, zpt=fits.getval(image_file, 'ZPTMAG'))
segmap = fits.getdata(segmap_file)

# Write out the sextable


for gf in glade_files:
    good_ids, glade_ids, glade_bmags = np.array([]), np.array([]), np.array([])

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

    # Output unique list of good sexcat IDs, glade ids, and glade bmags
    write_good_sexcat_ids(gf, image_file, good_ids, glade_ids, glade_bmags)

    # print(good_ids)
    # print(glade_ids)
    # print(glade_bmags)


t2 = time.time()
print("\n********************")
print("Execution time: %s" % (t2 - t1))
print("********************\n")