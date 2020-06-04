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


def write_good_sexcat_ids(glade_file, image_file, good_ids, glade_ids, glade_bmags, filtr, sex_mags, pixels, gal_coords):

    rows = []

    ascii_ecsv_fname = "%s_sexcat_good.txt" % glade_file.replace('.txt', '')
    ascii_ecsv_fpath = "%s/%s" % ("/data/LCO/Swope/logstch/gw190425/1/galaxy_demographics", ascii_ecsv_fname)
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
    region_fpath = "%s/%s.reg" % ("/data/LCO/Swope/logstch/gw190425/1/galaxy_demographics",
                                  glade_file.replace('.txt', ''))
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

    raise Exception("Stop!")

t1 = time.time()



# get all the swope files...
swope_files = []
swope_file_base_path = "/data/LCO/Swope/workstch/gw190425/1"
with open("./all_tiles_ascii.txt", 'r') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',')
    for row in csvreader:
        f = "%s/%s" % (swope_file_base_path, row[0])
        swope_files.append(f)


psf_shape = 31

for sf in swope_files:

    tokens = sf.split("/")[-1].split(".")
    field_name = tokens[0]
    photpipe_id = tokens[3].replace('_stch_1', '')


    glade_files = glob.glob('SwopeTiles/*%s*txt' % field_name)
    db_id = -9999
    for gf in glade_files:
        glade = at.Table.read(gf, format='ascii.ecsv')
        comment = glade.meta['comment']

        dcmp_photpipe_id = comment[0].split("=")[1].split(".")[3].replace("_stch_1", "")
        # import pdb; pdb.set_trace()

        if photpipe_id == dcmp_photpipe_id:
            db_id = comment[1].split("=")[1]
            break

    if db_id == -9999:
        raise Exception("Can't find match between glade files and `%s`!" % sf)


    glade_path = "./SwopeTiles"
    glade_file_name = "%s_%s.txt" % (db_id, field_name)
    glade_file_path = "%s/%s " % (glade_path, glade_file_name)

    dcmp_file = sf.replace('.fits', '.dcmp')
    mask_file = sf.replace('.fits', '.mask.fits.gz')



    # check if both files are on-disk
    # print(sf)
    # print(glade_file_path)
    # print(dcmp_file)
    # print(mask_file)

    if os.path.exists(sf) and os.path.exists(glade_file_path) \
            and os.path.exists(dcmp_file) and os.path.exists(mask_file):

        dcmp = txtobj(dcmp_file, cmpheader=True)
        dcmp_header = fits.getheader(dcmp_file)
        zpt = dcmp_header['ZPTMAG']
        filtr = dcmp_header['FILTER']

        mask_hdu = fits.open(mask_file)
        mask_data = mask_hdu[0].data.astype('float')

        # RUN SEXTRACTOR - this outputs the file *.sexcat
        segmap_file = sf.replace('.fits', '.check.fits')
        sextable = runsex(sf, segmapname=segmap_file, zpt=fits.getval(sf, 'ZPTMAG'))
        segmap = fits.getdata(segmap_file)

        good_ids, glade_ids, glade_bmags, measured_mags, gal_coords = [], [], [], [], []
        try:
            glade = at.Table.read(glade_file_path, format='ascii.ecsv')
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
            num_pix = segmap[good_galaxy_indices]
            pixels.append(len(num_pix))

        write_good_sexcat_ids(glade_file_name, sf, good_ids, glade_ids, glade_bmags, filtr, measured_mags, pixels,
                              gal_coords)

t2 = time.time()
print("\n********************")
print("Execution time: %s" % (t2 - t1))
print("********************\n")