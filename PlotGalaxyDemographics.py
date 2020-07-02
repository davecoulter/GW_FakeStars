import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import glob
import numpy as np
from astropy.io import fits
import os
from txtobj import txtobj
import re
from runsex import runsex
import sys
import astropy.table as at
# from astLib import astCoords
import astropy.coordinates as coord
from astropy import units as u
import time
from astropy.table import Table
import csv
import logging
from collections import OrderedDict


file_handler = logging.FileHandler(filename='PlotGalaxyDemographics.log')
stdout_handler = logging.StreamHandler(sys.stdout)
handlers = [file_handler, stdout_handler]
logging.basicConfig(
    level=logging.DEBUG,
    format='[%(asctime)s] {%(filename)s:%(lineno)d} %(levelname)s - %(message)s',
    handlers=handlers
)
logger = logging.getLogger('gw_fakes')
global_t1 = time.time()



measured_gal_bins = np.linspace(12.0, 19, 15)
pix_num_by_gal_mag = OrderedDict()
file_name_by_gal_mag = OrderedDict()

pix_num_by_file_by_mag = OrderedDict()
for g1, g2 in zip(measured_gal_bins[:-1], measured_gal_bins[1:]):
    pix_num_by_gal_mag[(g1, g2)] = 0.0
    file_name_by_gal_mag[(g1, g2)] = []
# print(pix_num_by_gal_mag.keys())

def get_gal_dict_key(gal_mag):
    return_key = None

    for k in pix_num_by_gal_mag.keys():
        bright = k[0]
        dim = k[1]
        if dim > gal_mag >= bright:
            return_key = k
            break

    if not return_key:
        logger.debug("Could not find key for gal_mag: %s" % gal_mag)
        raise Exception("Could not find key for gal_mag: %s" % gal_mag)

    return return_key


def generate_psf(psf_x, psf_xy, psf_y, psf_size):

    psf_model = np.zeros([psf_size, psf_size])
    for i in range(psf_size):
        for j in range(psf_size):
            x = i - psf_size/2
            y = j - psf_size/2
            zsq = 1/2.*(x**2./psf_x**2. + 2*psf_xy*x*y + y**2./psf_y**2.)
            psf_model[j, i] = (1 + zsq + 1/2.*zsq**2. + 1/6.*zsq**3.)**(-1)

    return psf_model


# Start with grabbing all SwopeDemographics files...
sexcat_files = glob.glob("./Fakes/Swope/SwopeDemographics/*_sexcat_good.txt")
# all_mags = []
all_pix_nums = []

for s in sexcat_files:
    sexcat_table = at.Table.read(s, format='ascii.ecsv')
    model_props = sexcat_table.meta['comment']
    f = model_props[0].split("=")[1].strip()

    for mag, num_pix in zip(list(sexcat_table["sex_mag"]), list(sexcat_table["num_pixels"])):
        key = get_gal_dict_key(mag)
        pix_num_by_gal_mag[key] += num_pix

        if f not in file_name_by_gal_mag[key]:
            file_name_by_gal_mag[key].append(f)

        if f not in pix_num_by_file_by_mag:
            pix_num_by_file_by_mag[f] = OrderedDict()

        if key not in pix_num_by_file_by_mag[f]:
            pix_num_by_file_by_mag[f][key] = 0.0

        pix_num_by_file_by_mag[f][key] += num_pix


# CREATE FILES BY GAL MAG -- as well as the accompanying template files...
# load template file
template_files = {}
with open("./Fakes/Swope/all_temps.txt", 'r') as tmp_file:
    csvreader = csv.reader(tmp_file, delimiter=' ', skipinitialspace=True)
    for row in csvreader:
        file_name = row[0]
        field_name = file_name.split(".")[0]
        template_files[field_name] = file_name

stars_needed = 3000
with open("./Fakes/Swope/Galaxy_Fakes/files_by_gal_mag.txt", 'w') as synopsis_file:
    for gal_key, file_list in file_name_by_gal_mag.items():

        total_fake_stars_per_gal_mag_bin = 0.0
        with open("./Fakes/Swope/Galaxy_Fakes/%s_%s_images.txt" % gal_key, 'w') as img_file, \
                open("./Fakes/Swope/Galaxy_Fakes/%s_%s_temps.txt" % gal_key, 'w') as tmp_file:
            for f in file_list:

                f_dcmp = f.split('/')[-1].replace('fits', 'dcmp')
                f_dcmp_header = fits.getheader("./Fakes/Swope/SwopeDemographics/" + f_dcmp)
                f_fwhm = f_dcmp_header['FWHM']
                f_pix_scale_arc_sec = f_dcmp_header['CD2_2'] * 3600.0
                f_psf_shape = int(np.ceil(2.0 * f_fwhm))

                file_name = f.split("/")[-1]
                img_file.write("%s\n" % file_name)

                field_name = file_name.split(".")[0]

                if field_name not in template_files:
                    print(f)
                    print("No matching template for field: %s" % field_name)
                    continue

                template = template_files[field_name]
                tmp_file.write("%s\n" % template)

                f_num_fake_stars = pix_num_by_file_by_mag[f][gal_key]/(f_psf_shape**2)
                total_fake_stars_per_gal_mag_bin += f_num_fake_stars

        num_needed = np.round(stars_needed / total_fake_stars_per_gal_mag_bin)
        if num_needed == 0:
            num_needed += 1

        synopsis_file.write("\nGalaxy Mag Bin (%0.1f, %0.1f); Fake Star Capacity: %0.2f; Runs: %0.2f\n" %
                      (gal_key[0], gal_key[1], total_fake_stars_per_gal_mag_bin, num_needed))



        for f in file_list:
            synopsis_file.write("%s\n" % f)
        synopsis_file.write("\n")

# region histogram plot
fig = plt.figure(figsize=(8, 8), dpi=300)
ax1 = fig.add_subplot(111)

hist_data = []
for g1, g2 in zip(measured_gal_bins[:-1], measured_gal_bins[1:]):
    bin_center = (g2-g1)/2.0 + g1
    num_pix = pix_num_by_gal_mag[(g1, g2)] / 961 # number of pixels in 1 fake star
    for n in range(int(num_pix)):
        hist_data.append(bin_center)
ax1.hist(hist_data, bins=measured_gal_bins)
ax1.set_xlabel('Galaxy Mag', fontsize=18)
# ax1.set_ylabel('Available Pixels', fontsize=18)
ax1.set_ylabel('Potential Fake Stars', fontsize=18)
fig.savefig('./Fakes/Swope/Galaxy_Fakes/SwopeGalaxyDemographics.png', bbox_inches='tight')
# endregion



# TESTS
do_tests = False
if do_tests:
    test_img = "s005aae0078.i.ut190425.1083_stch_1.sw.fits"
    fake_img = "s005aae0078.i.ut190425_fake.1083_stch_1.sw.fits"
    test_dcmp = "s005aae0078.i.ut190425.1083_stch_1.sw.dcmp"
    test_segmap = "s005aae0078.i.ut190425.1083_stch_1.sw.check.fits"
    test_mask = "s005aae0078.i.ut190425.1083_stch_1.sw.mask.fits.gz"

    image_hdu = fits.open(test_img)
    image_data = image_hdu[0].data.astype('float')
    segmap = fits.getdata(test_segmap)

    dcmp_header = fits.getheader(test_dcmp)
    zpt = dcmp_header['ZPTMAG']
    fwhm = dcmp_header['FWHM']

    pix_scale_arc_sec = dcmp_header['CD2_2']*3600.0
    psf_shape = int(np.ceil(2.0 * fwhm))
    print(psf_shape)

    mask_hdu = fits.open(test_mask)
    mask_data = mask_hdu[0].data.astype('float')

    # LOAD THE SERIALIZED PIXELS
    pix_by_sexcat_id = {}
    for sf in sexcat_files:
        glade = at.Table.read(sf, format='ascii.ecsv')
        comment = glade.meta['comment']
        fname = comment[0].split("=")[-1].split("/")[-1]

        if fname == test_img:

            # get pixel_good
            tokens = sf.split("_")

            pixel_good_file = "{}_{}_pixel_good.txt".format(tokens[0], tokens[1])
            pixel_good_table = at.Table.read(pixel_good_file, format='ascii.ecsv')

            sexcat_ids = list(pixel_good_table["sexcat_id"])
            good_pix_x = list(pixel_good_table["x"])
            good_pix_y = list(pixel_good_table["y"])

            for s_id in pixel_good_table["sexcat_id"]:
                if s_id not in pix_by_sexcat_id:
                    pix_by_sexcat_id[s_id] = [[], []]

            for i, s_id in enumerate(sexcat_ids):
                pix_by_sexcat_id[s_id][0].append(good_pix_x[i])
                pix_by_sexcat_id[s_id][1].append(good_pix_y[i])

            break

    # let's test by creating region boxes...
    with open("test_region_boxes.reg", 'w') as csvfile:

        csvfile.write("# Region file format: DS9 version 4.0 global\n\n")
        csvfile.write("global color=red\n")
        csvfile.write("image\n")

        for s_id, good_pix in pix_by_sexcat_id.items():
            all_x = good_pix[0]
            all_y = good_pix[1]

            width = np.max(all_x) - np.min(all_x)
            height = np.max(all_y) - np.min(all_y)

            central_x = np.min(all_x) + width/2.0
            central_y = np.min(all_y) + height/2.0

            csvfile.write('box(%s,%s,%s,%s,0) # width=1\n' % (central_x, central_y, width, height))

        print("Done w/ Region File")



    # Now test injecting point sources...
    psf_x, psf_xy, psf_y = dcmp_header['DPSIGX'],dcmp_header['DPSIGXY'], dcmp_header['DPSIGY']
    psf_model = generate_psf(psf_x, psf_xy, psf_y, psf_shape)

    test_mag = 20.8
    psf_mag = -2.5 * np.log10(np.sum(psf_model)) + zpt

    max_size = np.shape(psf_model)[0]
    dx = dy = int((max_size - 1) / 2)


    injected_fakes = []
    for s_id, good_pix in pix_by_sexcat_id.items():
        all_x = good_pix[0]
        all_y = good_pix[1]

        x_min = np.min(all_x)
        x_max = np.max(all_x)
        y_min = np.min(all_y)
        y_max = np.max(all_y)

        psf_loc_x = np.arange(x_min, x_max, psf_shape)
        psf_loc_y = np.arange(y_min, y_max, psf_shape)

        for x in psf_loc_x:
            for y in psf_loc_y:
                if segmap[y, x] == s_id:
                    injected_fakes.append((x, y, test_mag))

        for x, y, m in injected_fakes:
            psf_flux = 10 ** (-0.4 * (m - psf_mag))
            image_data[int(y) - dy:int(y) + dy + 1, int(x) - dx:int(x) + dx + 1] += psf_model * psf_flux

    image_hdu[0].data[:] = image_data
    image_hdu.writeto(fake_img, clobber=True, output_verify='ignore')

    # Write out Fakes
    fake_radius = (psf_shape/2.0)*pix_scale_arc_sec
    with open("fakes.reg", 'w') as csvfile:

        csvfile.write("# Region file format: DS9 version 4.0 global\n")
        csvfile.write("global color=green\n")
        csvfile.write("image\n")

        for x, y, m in injected_fakes:
            csvfile.write('circle(%s,%s,%s") # \n' % (x, y, fake_radius))
            csvfile.write('point(%s,%s) # point=cross\n' % (x, y))

    print("Done w/ Fakes Region File")

    global_t2 = time.time()
    print("\n********************")
    print("Execution time: %s" % (global_t2 - global_t1))
    print("********************\n")