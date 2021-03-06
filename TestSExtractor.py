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
import csv


def generate_psf(psf_x, psf_xy, psf_y, psf_size):
    psf_model = np.zeros([psf_size, psf_size])
    for i in range(psf_size):
        for j in range(psf_size):
            x = i - psf_size / 2
            y = j - psf_size / 2
            zsq = 1 / 2. * (x ** 2. / psf_x ** 2. + 2 * psf_xy * x * y + y ** 2. / psf_y ** 2.)
            psf_model[j, i] = (1 + zsq + 1 / 2. * zsq ** 2. + 1 / 6. * zsq ** 3.) ** (-1)

    return psf_model

def test_mask(arr_tup):

    arr_len = np.shape(arr_tup)[1]

    # Output region files as well
    region_fpath = "%s/%s.reg" % ("/data/LCO/Swope/logstch/gw190425/1/galaxy_demographics", "test")
    with open(region_fpath, 'w') as csvfile:
        csvfile.write("# Region file format: DS9 version 4.0 global\n\n")
        csvfile.write("global color=lightgreen\n")
        csvfile.write("image\n")

        for i in range(arr_len):

            x = arr_tup[1][i]
            y = arr_tup[0][i]
            csvfile.write('circle(%s,%s,1") # width=1\n' % (x, y))

        print("Done w/ Region File")

def write_good_sexcat_ids(glade_file, image_file, good_ids, glade_ids, glade_bmags, filtr, sex_mags, pixels, gal_coords):

    rows = []

    ascii_ecsv_fname = "%s_sexcat_good.txt" % glade_file.replace('.txt', '')
    ascii_ecsv_fpath = "%s/%s" % (".", ascii_ecsv_fname)
    print("Creating `%s`" % ascii_ecsv_fpath)

    # Build ascii.ecsv formatted output
    cols = ['sexcat_id', 'glade_id', 'ra_dec', 'dec_dec', 'glade_B', 'filter', 'sex_mag', 'num_pixels']
    dtype = ['i4', 'i4', 'f8', 'f8', 'f8', 'U64', 'f8', 'i4']
    result_table = Table(dtype=dtype, names=cols)
    meta = ["{key}={value}".format(key="image_file", value=image_file)]
    result_table.meta['comment'] = meta

    for sexcat_id, glade_id, coord_tup, b, sex_mag, num_pix in \
            zip(good_ids, glade_ids, gal_coords, glade_bmags, sex_mags, pixels):
        rows.append([sexcat_id, glade_id, coord_tup[0], coord_tup[1], b, filtr, sex_mag, num_pix])

    for r in rows:
        result_table.add_row(r)

    result_table.write(ascii_ecsv_fpath, overwrite=True, format='ascii.ecsv')

    # Output region files as well
    region_fpath = "%s/%s.reg" % (".", glade_file.replace('.txt', ''))
    with open(region_fpath, 'w') as csvfile:

        csvfile.write("# Region file format: DS9 version 4.0 global\n\n")
        csvfile.write("global color=lightgreen\n")
        csvfile.write("ICRS\n")

        for r in rows:
            glade_id = r[1]
            ra = r[2]
            dec = r[3]
            csvfile.write('circle(%s,%s,30") # width=2 text="%s"\n' % (ra, dec, glade_id))

        print("Done w/ Region File")

t1 = time.time()

print("running SExtractor")

psf_shape = 31

dcmp_file = "s005aae0078.i.ut190425.1083_stch_1.sw.dcmp"
dcmp = txtobj(dcmp_file, cmpheader=True)
dcmp_header = fits.getheader(dcmp_file)
zpt = dcmp_header['ZPTMAG']
filtr = dcmp_header['FILTER']

mask_file = "s005aae0078.i.ut190425.1083_stch_1.sw.mask.fits.gz"
mask_hdu = fits.open(mask_file)
mask_data = mask_hdu[0].data.astype('float')

image_file = "s005aae0078.i.ut190425.1083_stch_1.sw.fits"
field_name = image_file.split('.')[0]
glade_files = glob.glob('*%s*txt' % field_name)
gf = glade_files[0]

# RUN SEXTRACTOR - this outputs the file *.sexcat
segmap_file = image_file.replace('.fits', '.check.fits')
sextable = runsex(image_file, segmapname=segmap_file, zpt=fits.getval(image_file, 'ZPTMAG'))
segmap = fits.getdata(segmap_file)


good_ids, glade_ids, glade_bmags, measured_mags, gal_coords = [], [], [], [], []
try:
    glade = at.Table.read(gf, format='ascii.ecsv')
except:
    raise Exception("Can't read file...")

# Identify matches between teglon/GLADE and what source extractor detects
# keep a list of matches ("good_ids")
for i, sex_num in enumerate(sextable.NUMBER):

    # separation between sex source and all glade galaxies in arc sec
    sep = astCoords.calcAngSepDeg(glade['Galaxy_RA'],
                                  glade['Galaxy_Dec'],
                                  sextable.X_WORLD[i],
                                  sextable.Y_WORLD[i])*3600.



    # sanity
    if len(sep) == 0:
        continue

    min_sep = np.min(sep)
    closest_sep_mask = sep == min_sep

    # if within 2 arc sec. # No duplicates
    if min_sep <= 2 and sex_num not in good_ids:

        # import pdb; pdb.set_trace()

        good_ids.append(sex_num)
        glade_ids.append(glade['Galaxy_ID'][closest_sep_mask][0])
        glade_bmags.append(glade['B'][closest_sep_mask][0])
        measured_mags.append(sextable.MAG_AUTO[i])
        gal_coords.append((sextable.X_WORLD[i], sextable.Y_WORLD[i]))

# Zero out entries in the segmap that not matched galaxies
for i in sextable.NUMBER:
    if i not in good_ids:
        segmap[segmap == i] = 0

# get the indices (# of pixels) for each galaxy based on sextractor
pixels = []
for i in good_ids:
    good_galaxy_indices = np.where((mask_data != 144.0) & (segmap == i))

    test_mask(good_galaxy_indices)
    raise Exception("Stop")

    num_pix = segmap[good_galaxy_indices]
    pixels.append(len(num_pix))

# psf_x, psf_xy, psf_y = dcmp_header['DPSIGX'],dcmp_header['DPSIGXY'], dcmp_header['DPSIGY']
# psf_model = generate_psf(psf_x, psf_xy, psf_y, psf_shape)

write_good_sexcat_ids(gf, image_file, good_ids, glade_ids, glade_bmags, filtr, measured_mags, pixels, gal_coords)

    # print(good_ids)
    # print(glade_ids)
    # print(glade_bmags)


t2 = time.time()
print("\n********************")
print("Execution time: %s" % (t2 - t1))
print("********************\n")