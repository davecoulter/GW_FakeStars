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
import logging
from astropy.modeling.models import Sersic1D, Sersic2D
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors

logging.basicConfig(filename='GalaxyDemographics.log', level=logging.DEBUG)


# File paths
# # SWOPE
# file_base_path = "/data/LCO/Swope/workstch/gw190425/1"
# all_tiles_asci = "./Fakes/Swope/all_tiles_ascii.txt"
# template_file_base_path = "/data/LCO/Swope/workstch/gw190425tmpl/1"
# all_temps = "./Fakes/Swope/all_temps.txt"
# tile_file_pattern_match = "./Fakes/Swope/SwopeTiles/*%s*txt"
# glade_path = "./Fakes/Swope/SwopeTiles"
# demographics_path = "/data/LCO/Swope/logstch/gw190425/1/galaxy_demographics"

# THACHER
file_base_path = "/data2/THACHER/workspace/gw190425/1"
all_tiles_asci = "./Fakes/Thacher/all_tiles_ascii.txt"
template_file_base_path = "/data2/THACHER/workspace/gw190425tmpl/1"
all_temps = "./Fakes/Thacher/all_temps.txt"
tile_file_pattern_match = "./Fakes/Thacher/ThacherTiles/*%s*txt"
glade_path = "./Fakes/Thacher/ThacherTiles"
demographics_path = "/data2/THACHER/logs/gw190425/1/galaxy_demographics"


def write_good_sexcat_ids(glade_file, image_file, good_ids, glade_ids, glade_bmags, filtr, sex_mags, pixels, gal_coords,
                          pixel_tuple_dict, flux_radii, gal_xy, cxx, cyy, cxy, A, B, theta, ellip):

    galaxy_rows = []

    # Tile synoptic information
    ascii_ecsv_fname = "%s_sexcat_good.txt" % glade_file.replace('.txt', '')
    ascii_ecsv_fpath = "%s/%s" % (demographics_path, ascii_ecsv_fname)
    print("Creating `%s`" % ascii_ecsv_fpath)
    cols = ['sexcat_id', 'glade_id', 'ra_dec', 'dec_dec', 'glade_B', 'filter', 'sex_mag', 'num_pixels',
            'flux_radius', 'x', 'y', 'cxx', 'cyy', 'cxy', 'A', 'B', 'theta', 'ellip']
    dtype = ['i4', 'i4', 'f8', 'f8', 'f8', 'U64', 'f8', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8']
    result_table = Table(dtype=dtype, names=cols)
    meta = ["{key}={value}".format(key="image_file", value=image_file)]
    result_table.meta['comment'] = meta

    for sexcat_id, glade_id, coord_tup, b, sex_mag, num_pix, flux_rad, xy_tup, _cxx, _cyy, _cxy, _A, _B, _theta, _ellip in \
            zip(good_ids, glade_ids, gal_coords, glade_bmags, sex_mags, pixels, flux_radii, gal_xy, cxx, cyy, cxy, A, B, theta, ellip):
        galaxy_rows.append([sexcat_id, glade_id, coord_tup[0], coord_tup[1], b, filtr, sex_mag, num_pix, flux_rad, xy_tup[0], xy_tup[1], _cxx, _cyy, _cxy, _A, _B, _theta, _ellip])

    for r in galaxy_rows:
        result_table.add_row(r)

    result_table.write(ascii_ecsv_fpath, overwrite=True, format='ascii.ecsv')


    # Tile good pix
    ascii_ecsv_fname2 = "%s_pixel_good.txt" % glade_file.replace('.txt', '')
    ascii_ecsv_fpath2 = "%s/%s" % (demographics_path, ascii_ecsv_fname2)
    print("Creating `%s`" % ascii_ecsv_fpath2)
    # Build ascii.ecsv formatted output
    cols2 = ['sexcat_id', 'x', 'y', 'sep', 'weight']
    dtype2 = ['i4', 'i4', 'i4', 'f8', 'f8']
    result_table2 = Table(dtype=dtype2, names=cols2)

    weighted_pixels = {}

    print("Num galaxies: %s" % len(pixel_tuple_dict))
    for sxct_id, pixel_tuple in pixel_tuple_dict.items():

        # import pdb; pdb.set_trace()
        weighted_pixels[sxct_id] = []
        # Build table for valid galaxy pixels

        # pixel_tuple contains:
        #   [0] good_galaxy_indices
        #   [1] galaxy X/Y tuple
        #   [2] galaxy flux radius
        ggi = pixel_tuple[0]
        gxy = pixel_tuple[1]
        gfr = pixel_tuple[2]
        gel = pixel_tuple[3]
        gth = pixel_tuple[4]

        print("sxct_id: %s; num_pix: %s" % (sxct_id, len(ggi[0])))

        # assuming this... from Ryan (7/15/2020): "it looks like n = 2 might be a good sersic index it is both
        # between spirals and ellipticals and it avoids some issues with smaller galaxies that are only a few times
        # the size of your PSF"
        # https: // sextractor.readthedocs.io / en / latest / Position.html  # fig-ellipse
        # http: // astroa.physics.metu.edu.tr / MANUALS / sextractor / Guide2source_extractor.pdf
        s2D = Sersic2D(amplitude=1, r_eff=gfr, n=2.0, x_0=gxy[0], y_0=gxy[1], ellip=gel, theta=np.radians(gth))

        arr_len = np.shape(ggi)[1]
        for i in range(arr_len):
            x = ggi[1][i]
            y = ggi[0][i]

            sep = np.sqrt((x - gxy[0])**2.0 + (y - gxy[1])**2.0)
            # weight = s1D(sep)
            weight = s2D(x, y)

            weighted_pixels[sxct_id].append((x, y, sep, weight))
            result_table2.add_row([sxct_id, x, y, sep, weight])
    print("Writing good pixels")
    result_table2.write(ascii_ecsv_fpath2, overwrite=True, format='ascii.ecsv')

    # Tile region
    region_fpath = "%s/%s.reg" % (demographics_path,
                                  glade_file.replace('.txt', '_galaxies'))
    with open(region_fpath, 'w') as csvfile:

        csvfile.write("# Region file format: DS9 version 4.0 global\n\n")
        csvfile.write("global color=lightgreen\n")
        csvfile.write("ICRS\n")

        for r in galaxy_rows:
            sexcat_id = r[0]
            glade_id = r[1]
            ra = r[2]
            dec = r[3]
            csvfile.write('circle(%s,%s,30") # width=2 text="sxct_id: %s | db_id: %s"\n' % (ra, dec, int(sexcat_id), int(glade_id)))

        print("Done w/ Galaxy Position Region File")


    region_fpath2 = "%s/%s.reg" % (demographics_path,
                                  glade_file.replace('.txt', '_weighted_pixels'))
    with open(region_fpath2, 'w') as csvfile:

        csvfile.write("# Region file format: DS9 version 4.0 global\n\n")
        csvfile.write("global color=red\n")
        csvfile.write("image\n")

        for sxct_id, pixel_list in weighted_pixels.items():

            weights = []
            for pixel_tuple in pixel_list:
                weights.append(pixel_tuple[3])

            if len(weights) == 0:
                print("No weights for sexcat_id: %s" % sxct_id)
                logging.debug("No weights for sexcat_id: %s" % sxct_id)
                continue

            norm = colors.LogNorm(min(weights), max(weights))

            for pixel_tuple in pixel_list:
                # (x, y, sep, weight)
                x = pixel_tuple[0]
                y = pixel_tuple[1]
                w = pixel_tuple[3]
                clr = plt.cm.viridis(norm(w))
                hex_clr = matplotlib.colors.rgb2hex(clr)

                csvfile.write('box(%s,%s,%s,%s) # width=2 color=%s\n' % (x, y, 1.0, 1.0, hex_clr))

        print("Done w/ Pixel Weight Region File")

def test_mask(arr_tup):

    arr_len = np.shape(arr_tup)[1]

    # Output region files as well
    region_fpath = "%s/%s.reg" % (demographics_path, "test")
    with open(region_fpath, 'w') as csvfile:
        csvfile.write("# Region file format: DS9 version 4.0 global\n\n")
        csvfile.write("global color=lightgreen\n")
        csvfile.write("image\n")

        for i in range(arr_len):

            x = arr_tup[1][i]
            y = arr_tup[0][i]
            csvfile.write('circle(%s,%s,1) # width=1\n' % (x, y))

        print("Done w/ Region File")

global_t1 = time.time()
# get all the swope files...
img_files = []
with open(all_tiles_asci, 'r') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',')
    for row in csvreader:
        f = "%s/%s" % (file_base_path, row[0])
        img_files.append(f)

template_files = []
with open(all_temps, 'r') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',')
    for row in csvreader:
        f = "%s/%s" % (template_file_base_path, row[0])
        template_files.append(f)


num_sf = len(img_files)
for img_file_index, img_file in enumerate(img_files):

    print("\n********************")
    print("Processing %s/%s..." % (img_file_index + 1, num_sf))
    print("********************\n")
    t1 = time.time()

    tokens = img_file.split("/")[-1].split(".")
    field_name = tokens[0]
    photpipe_id = tokens[3].replace('_stch_1', '')

    glade_files = glob.glob(tile_file_pattern_match % field_name)
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
        print("Can't find match between glade files and `%s`!" % img_file)
        print("Skipping %s" % img_file)
        logging.debug("Can't find match between glade files and `%s`!" % img_file)
        logging.debug("Skipping %s" % img_file)
        continue

    glade_file_name = "%s_%s.txt" % (db_id, field_name)
    glade_file_path = "%s/%s" % (glade_path, glade_file_name)

    # check if this file has already been processed
    ascii_ecsv_fname = "%s_sexcat_good.txt" % glade_file_name.replace('.txt', '')
    ascii_ecsv_fpath = "%s/%s" % (demographics_path, ascii_ecsv_fname)
    if os.path.exists(ascii_ecsv_fpath):
        print("%s already processed! skipping..." % glade_file_name)
        logging.debug("%s already processed! skipping..." % glade_file_name)
        continue



    img_dcmp_file = img_file.replace('.fits', '.dcmp')
    img_mask_file = img_file.replace('.fits', '.mask.fits.gz')

    # get template file path from field name
    tmp_file = None
    for tf in template_files:
        if field_name in tf:
            tmp_file = tf
            break

    if tmp_file is None:
        raise Exception("No matching template file for image: %s" % img_file)

    tmp_dcmp_file = tmp_file.replace('.fits', '.dcmp')
    tmp_mask_file = tmp_file.replace('.fits', '.mask.fits.gz')
    if not os.path.exists(tmp_mask_file):
        tmp_mask_file = tmp_file.replace('.fits', '.mask.fits')

    if os.path.exists(img_file) and os.path.exists(glade_file_path) \
            and os.path.exists(img_dcmp_file) and os.path.exists(img_mask_file) \
            and os.path.exists(tmp_file) and os.path.exists(tmp_dcmp_file) and os.path.exists(tmp_mask_file):

        img_dcmp = txtobj(img_dcmp_file, cmpheader=True)
        img_dcmp_header = fits.getheader(img_dcmp_file)
        img_zpt = img_dcmp_header['ZPTMAG']
        filtr = img_dcmp_header['FILTER']

        img_mask_hdu = fits.open(img_mask_file)
        img_mask_data = img_mask_hdu[0].data.astype('float')

        tmp_mask_hdu = fits.open(tmp_mask_file)
        tmp_mask_data = tmp_mask_hdu[0].data.astype('float')

        # RUN SEXTRACTOR - this outputs the file *.sexcat
        segmap_file = img_file.replace('.fits', '.check.fits')
        sextable = runsex(img_file, segmapname=segmap_file, zpt=fits.getval(img_file, 'ZPTMAG'))
        segmap = fits.getdata(segmap_file)

        good_ids, glade_ids, glade_bmags, measured_mags, gal_coords, gal_xy, flux_radii, cxx, cyy, cxy, A, B, theta, ellip = \
            [], [], [], [], [], [], [], [], [], [], [], [], [], []
        try:
            glade = at.Table.read(glade_file_path, format='ascii.ecsv')

            if len(glade['Galaxy_RA']) == 0:
                print("No glade information in this field, skipping!")
                logging.debug("No glade information in this field, skipping!")
                continue
        except:
            logging.debug("Can't read `%s`!" % glade_file_path)
            logging.debug("Exiting!")
            raise Exception("Can't read file `%s`" % glade_file_path)

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
                logging.debug("No sep for sextable.NUMBER=%s for `%s` and `%s`" % (sex_num, glade_file_path, img_file))
                logging.debug("Skipping sextable.NUMBER=%s" % sex_num)
                print("No sep for sextable.NUMBER=%s for `%s` and `%s`" % (sex_num, glade_file_path, img_file))
                print("Skipping sextable.NUMBER=%s" % sex_num)
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
                flux_radii.append(sextable.FLUX_RADIUS[i])
                gal_xy.append((sextable.X_IMAGE[i], sextable.Y_IMAGE[i]))
                cxx.append(sextable.CXX_IMAGE[i])
                cyy.append(sextable.CYY_IMAGE[i])
                cxy.append(sextable.CXY_IMAGE[i])
                A.append(sextable.A_IMAGE[i])
                B.append(sextable.B_IMAGE[i])
                theta.append(sextable.THETA_IMAGE[i])
                ellip.append(sextable.ELLIPTICITY[i])

                # Zero out entries in the segmap that not matched galaxies
        for i in sextable.NUMBER:
            if i not in good_ids:
                segmap[segmap == i] = 0

        # get the indices (# of pixels) for each galaxy based on sextractor
        pixels = []
        pixel_tup_dict = {}
        for i, good_id in enumerate(good_ids):

            # CHANGE THIS TO ALSO RESPECT THE TEMPLATE MASK!
            good_galaxy_indices = np.where((img_mask_data != 144.0) & (tmp_mask_data != 144.0) & (segmap == good_id))

            # print("sxct_id: %s" % good_id)
            # print("\tnum good indices: %s" % len(good_galaxy_indices[0]))
            # import pdb; pdb.set_trace()

            # send over the good pixel indices, the galaxy X/Y position, and the galaxy flux radius
            pixel_tup_dict[good_id] = (good_galaxy_indices, gal_xy[i], flux_radii[i], ellip[i], theta[i])
            gal_pix = segmap[good_galaxy_indices]
            pixels.append(len(gal_pix))

        write_good_sexcat_ids(glade_file_name, img_file, good_ids, glade_ids, glade_bmags, filtr, measured_mags, pixels,
                              gal_coords, pixel_tup_dict, flux_radii, gal_xy, cxx, cyy, cxy, A, B, theta, ellip)

        t2 = time.time()
        print("\n********************")
        print("Done processing %s/%s: %s" % (img_file_index + 1, num_sf, (t2 - t1)))
        print("********************\n")
    else:

        if not os.path.exists(img_file):
            print("Path doesn't exist for: `%s`" % img_file)
            logging.debug("Path doesn't exist for: `%s`" % img_file)
        if not os.path.exists(glade_file_path):
            print("Path doesn't exist for: `%s`" % glade_file_path)
            logging.debug("Path doesn't exist for: `%s`" % glade_file_path)
        if not os.path.exists(img_dcmp_file):
            print("Path doesn't exist for: `%s`" % img_dcmp_file)
            logging.debug("Path doesn't exist for: `%s`" % img_dcmp_file)
        if not os.path.exists(img_mask_file):
            print("Path doesn't exist for: `%s`" % img_mask_file)
            logging.debug("Path doesn't exist for: `%s`" % img_mask_file)

        if not os.path.exists(tmp_file):
            print("Path doesn't exist for: `%s`" % tmp_file)
            logging.debug("Path doesn't exist for: `%s`" % tmp_file)
        if not os.path.exists(tmp_dcmp_file):
            print("Path doesn't exist for: `%s`" % tmp_dcmp_file)
            logging.debug("Path doesn't exist for: `%s`" % tmp_dcmp_file)
        if not os.path.exists(tmp_mask_file):
            print("Path doesn't exist for: `%s`" % tmp_mask_file)
            logging.debug("Path doesn't exist for: `%s`" % tmp_mask_file)

        print("Skipping %s" % img_file)
        logging.debug("Skipping %s" % img_file)

    # if sf_index == 0:
    #     print("\n\nDebug stop!!\n\n")
    #     break


global_t2 = time.time()
print("\n********************")
print("Execution time: %s" % (global_t2 - global_t1))
print("********************\n")