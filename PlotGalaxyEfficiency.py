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
from Binomial import *
from scipy.special import erf
from scipy.optimize import minimize


# GLOBAL SETTINGS
#   These settings control global constants/formatters
base_path = "./Fakes/Swope/Galaxy_Fakes"
server_img_path = "/data/LCO/Swope/workstch/gw190425/1"

#   These settings control the x-axis resolution of the efficiency curve
bin_size = 0.2 # mag
bright = 18
dim = 25
efficiency_magbins = np.linspace(18, 25, (dim-bright)/bin_size)

# TODO: this process needs to be formalized to see what is what in band, and only process the fakes within 1 band...
current_band = "i"

# RUN SETTINGS
#   Which galaxy mag bin to process, and which data files to use

gal_paths = [
    "13.0_13.5",
    "13.5_14.0",
    "14.0_14.5",
    "14.5_15.0",
    "15.0_15.5",
    "15.5_16.0",
    "16.0_16.5",
    "16.5_17.0",
    "17.0_17.5"
]
gal_subdirs = [
    [1,2,3],
    [4,5,6],
    [7,8,9],
    [10,11,12],
    [13,14,15],
    [16,17,18],
    [19,20,21],
    [22,23,24],
    [25,26,27]
]

dcmp_type = "cut.dcmp"
# dcmp_type = "diff.dcmp"

snr = 10.0
model_mags = np.linspace(18.0, 25.0, 500)
# unrecovered_range = (18, 25)
unrecovered_range = (19, 20)


class GalEfficiency:

    measured_fakes_file = "measuredfakemags.txt"
    efficiency_file = "efficiencies_{}.txt"  # formatted with SNR

    def __init__(self, base_gal_bin_path, efficiency_magbins, snr, dcmp_type):

        self.dcmp_type = dcmp_type
        self.base_gal_bin_path = base_gal_bin_path
        self.output_dir = "{}/{}.{}".format(base_path, base_gal_bin_path, dcmp_type)
        # Create output directory if it does not exist
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.efficiency_magbins = efficiency_magbins

        # Initialize measurefakemags output file
        self.measured_fakes_output = "{}/{}".format(self.output_dir, GalEfficiency.measured_fakes_file)
        with open(self.measured_fakes_output, 'w') as fout:
            print('# dcmpfile x y sim_mag mag_gal mag_gal_GLADE GLADE_gal_id detmag detmagerr snr weight', file=fout)

        self.snr = snr
        self.efficiency_output = "{}/{}".format(self.output_dir, GalEfficiency.efficiency_file.format(snr))

        self.unweighted_efficiency = []
        self.weighted_efficiency = []
        self.number = []
        self.central_mag_bins = []
        self.error_high = []
        self.error_low = []

        # Initialized by call to compute_efficiency()
        self.mean_errors = None
        self.fitted_x0 = None
        self.fitted_a = None
        self.fitted_n = None

    def chi_square_s_curve(self, params, mag_bins=None, obs_efficiency=None, errors=None):
        x0, a, n = params
        model_efficiency = self.limiting_mag_model(mag_bins, x0, a, n)
        chi_square = np.sum(((model_efficiency - obs_efficiency) / errors) ** 2.0)
        return chi_square

    def compute_efficiency(self, initial_guess):

        measuredfakemags = txtobj(self.measured_fakes_output)

        # Write file out so we can inspect it for debug purposes
        with open(self.efficiency_output, 'w') as fout:
            # print('# mag deteff N mean_gal_mag', file=fout)
            print('# mag unweighted_eff weighted_eff N mean_gal_mag', file=fout)

            mag_sims = measuredfakemags.sim_mag
            detmags = measuredfakemags.detmag
            mag_gal = measuredfakemags.mag_gal
            weights = measuredfakemags.weight
            snr = measuredfakemags.snr

            for m1, m2 in zip(self.efficiency_magbins[:-1], self.efficiency_magbins[1:]):

                iDet = np.where((mag_sims[measuredfakemags.snr >= self.snr] >= m1) &
                                (mag_sims[measuredfakemags.snr >= self.snr] <= m2) &
                                (detmags[measuredfakemags.snr >= self.snr] != -99))[0]

                # iDet = np.where((mag_sims >= m1) & (mag_sims <= m2) & (detmags != -99) &
                #                 (snr >= self.snr))[0]

                detWeight = weights[iDet]

                iAll = np.where((mag_sims >= m1) & (mag_sims <= m2))[0]
                allWeight = weights[iAll]

                iAll2 = np.where((mag_sims >= m1) & (mag_sims <= m2) &
                                 (detmags == -99))[0]

                mean_gal_mag = np.mean(mag_gal[iDet], axis=0)
                undetected_mean_gal_mag = np.mean(mag_gal[iAll2], axis=0)

                if len(iAll):

                    # weighted_eff = np.sum(detWeight)/np.sum(allWeight)
                    unweighted_eff = len(iDet) / float(len(iAll))
                    weighted_eff = len(iDet) / float(len(iAll)) # test

                    # print('%.2f %.3f %0.3f %s %0.3f' % ((m1 + m2) / 2., unweighted_eff, weighted_eff, len(iAll), mean_gal_mag), file=fout)
                    if len(iDet) > 0:
                        # print('%.2f %.3f %s %.3f %.3f' % ((m1 + m2) / 2., len(iDet) / float(len(iAll)), len(iDet) / float(len(iAll)), len(iAll), mean_gal_mag),
                        #       file=fout)
                        print('%.2f %.3f %0.3f %s %0.3f' % (
                        (m1 + m2) / 2., unweighted_eff, weighted_eff, len(iAll), mean_gal_mag), file=fout)

                    else:
                        # print('%.2f %.3f %s %.3f %.3f' % ((m1 + m2) / 2., len(iDet) / float(len(iAll)), len(iDet) / float(len(iAll)), len(iAll), undetected_mean_gal_mag),
                        #       file=fout)
                        print('%.2f %.3f %.3f %s %.3f' % (
                        (m1 + m2) / 2., unweighted_eff, weighted_eff, len(iAll), undetected_mean_gal_mag), file=fout)
                else:
                    print('%.2f %.3f %0.3f %s %.3f' % ((m1 + m2) / 2., np.nan, np.nan, len(iAll), np.nan), file=fout)

        # Read in file. In the future, we can just go directly here without the writing step...
        with open(self.efficiency_output, 'r') as fin:
            csv_reader = csv.reader(fin, delimiter=' ')
            next(csv_reader)

            for row in csv_reader:
                # if float(row[2]) > 0:
                number = float(row[3])
                if number > 0:
                    self.central_mag_bins.append(float(row[0]))
                    self.unweighted_efficiency.append(float(row[1]))
                    self.weighted_efficiency.append(float(row[2]))
                    self.number.append(number)

        for eff, num in zip(self.unweighted_efficiency, self.number):
            yes = eff * num
            no = (1.0 - eff) * num
            ci = binomial(yes, no)

            self.error_high.append(ci[0])
            self.error_low.append(ci[1])

        self.mean_errors = np.mean([self.error_high, self.error_low], axis=0)
        minimize_result = minimize(self.chi_square_s_curve, initial_guess,
                                   args=(self.central_mag_bins, self.weighted_efficiency, self.mean_errors))

        self.fitted_x0 = minimize_result.x[0]
        self.fitted_a = minimize_result.x[1]
        self.fitted_n = minimize_result.x[2]

        # Write out model parameter result
        efficiency_param_file = "{}/{}".format(gal_efficiency_3.output_dir, "efficiency_params.txt")
        with open(efficiency_param_file, 'w') as csvfile:
            csvfile.write("# x0 a n\n")
            csvfile.write("%s %s %s" % (self.fitted_x0, self.fitted_a, self.fitted_n))

    # Limiting mag model uses Eq 12 from:
    #   http://www.cfht.hawaii.edu/Science/CFHLS/T0007/T0007-docsu12.html
    #   I have removed the 100x prefactor to keep in decimal percent, and introduced a normalization `n` that allows
    #   for the maximum efficiency to be different than 1.0
    #   Parameters:
    #       x0 = turnover in function
    #       a = slope at turnover
    #       n = normalization
    def limiting_mag_model(self, input_mags, x0, a, n):

        y = n * (1.0 - (erf(a * (input_mags - x0)) + 1.0) / 2.0)
        return y

    def computed_lim_mag_model(self, input_mags):
        # Sanity
        if self.fitted_x0 is None or self.fitted_a is None or self.fitted_n is None:
            raise Exception("Model not initialized yet! Must invoke GalEfficiency.compute_efficiency()! Exiting...")

        return self.limiting_mag_model(input_mags, self.fitted_x0, self.fitted_a, self.fitted_n)

    ### Version that allows a shift in input mags...
    # def limiting_mag_model(input_mags, x0, a, n, d_mag=0.0):
    #     delta_mags = np.zeros(len(input_mags)) + d_mag
    #     y = n*(1.0 - (erf(a*((input_mags + delta_mags) - x0)) + 1.0)/2.0) # Removed the 100x prefactor
    #     return y


class GalFake:
    fakes_file = "fakemags.txt"

    def __init__(self, base_gal_bin_path, gal_bin_subdir):

        self.gal_bin_subdir = gal_bin_subdir
        self.fakes_path = "{}/{}_{}/{}".format(base_path,
                                               base_gal_bin_path,
                                               gal_bin_subdir,
                                               GalFake.fakes_file)
        self.dcmp_path = "{}/{}_{}".format(base_path,
                                           base_gal_bin_path,
                                           gal_bin_subdir)

    def measure_fakemags(self, fake_txtobj, dcmp_txtobj, reconstructed_filepath, zpt, output_file, weighted_pixel_dict,
                         unrecovered_reg_file=None, dcmp_file=None, unrecovered_range=None,
                         ds9_cmd_file=None):


        dcmp_header = fits.getheader(dcmp_file)
        gal_fake_fwhm_factor = 3.0
        fwhm = dcmp_header['FWHM']
        pix_scale_arc_sec = np.abs(dcmp_header['CD2_2']) * 3600.0
        psf_shape = int(np.ceil(gal_fake_fwhm_factor * fwhm))
        fake_radius = (psf_shape / 2.0) * pix_scale_arc_sec
        m1 = None
        m2 = None
        empty_reg = True
        if unrecovered_range:
            m1 = unrecovered_range[0]
            m2 = unrecovered_range[1]

        # initialize region file...
        with open(unrecovered_reg_file, 'w') as csvfile:
            csvfile.write("# Region file format: DS9 version 4.0 global\n")
            csvfile.write('global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
            csvfile.write("image\n")

        # print(reconstructed_filepath)
        for x, y, mag, mag_gal, mag_gal_GLADE, GLADE_gal_id \
                in zip(fake_txtobj.x[fake_txtobj.dcmpfile == reconstructed_filepath],
                       fake_txtobj.y[fake_txtobj.dcmpfile == reconstructed_filepath],
                       fake_txtobj.mag[fake_txtobj.dcmpfile == reconstructed_filepath],
                       fake_txtobj.mag_gal[fake_txtobj.dcmpfile == reconstructed_filepath],
                       fake_txtobj.mag_gal_GLADE[fake_txtobj.dcmpfile == reconstructed_filepath],
                       fake_txtobj.GLADE_gal_id[fake_txtobj.dcmpfile == reconstructed_filepath]):

            sep = np.sqrt((x - dcmp_txtobj.Xpos) ** 2. + (y - dcmp_txtobj.Ypos) ** 2.)
            w = weighted_pixel_dict[(x, y)]

            # find where fake mag injection is < 2 pixels from DCMP measured source
            # sep_arr = np.where(sep < 2)[0]
            # sep_arr = np.where(sep < 5)[0]
            if len(np.where(sep < 5)[0]) > 0:

                # source confusion: skip if there are multiple matches within 2 pix
                iMin = np.where(sep == np.min(sep))[0]
                if len(iMin) > 1:
                    continue

                # HACK: if there's only one value in dcmp_txtobj.flux wrap it in an array so indexing doesn't fail.
                flux_val = None
                dflux_val = None
                try:
                    flux_val = dcmp_txtobj.flux[iMin]
                    dflux_val = dcmp_txtobj.dflux[iMin]
                except:
                    flux_val = np.asarray([dcmp_txtobj.flux])[iMin]
                    dflux_val = np.asarray([dcmp_txtobj.dflux])[iMin]

                # Sanity on the values...
                if (flux_val < 0.0):  # reject if flux is < 0
                    with open(output_file, 'a') as fout:
                        print('%s %.2f %.2f %.3f %.3f %.3f %.0f %.0f %.0f %0.f %.3f' %
                              (reconstructed_filepath, x, y, mag, mag_gal, mag_gal_GLADE, GLADE_gal_id, -99, -99, -99,
                               w), file=fout)
                else:
                    with open(output_file, 'a') as fout:
                        print('%s %.2f %.2f %.3f %.3f %.3f %.0f %.3f %.3f %.3f %.3f' %
                              (reconstructed_filepath, x, y, mag, mag_gal, mag_gal_GLADE, GLADE_gal_id,
                               -2.5 * np.log10(flux_val) + zpt,  # calculate mag
                               1.086 * dflux_val / flux_val,  # calculate mag err
                               flux_val / dflux_val, # SNR
                               w), file=fout)
            else:
                with open(output_file, 'a') as fout:
                    print('%s %.2f %.2f %.3f %.3f %.3f %.0f %.0f %.0f %0.f %.3f' %
                          (reconstructed_filepath, x, y, mag, mag_gal, mag_gal_GLADE, GLADE_gal_id, -99, -99, -99, w),
                          file=fout)

                # If this star is within the mag range we're considering...
                if mag >= m1 and mag <= m2:
                    empty_reg = False
                    with open(unrecovered_reg_file, 'a') as csvfile:
                        csvfile.write('circle(%s,%s,%s") # color=blue width=4\n' % (x, y, fake_radius))
                        csvfile.write('# text(%s,%s) textangle=360 color=blue width=3 font="helvetica 10 bold roman" text={%0.1f} \n' % (x, y, mag))

        if empty_reg:
            # delete region file
            os.remove(unrecovered_reg_file)
        else:
            if ds9_cmd_file:
                with open(ds9_cmd_file, 'a') as csvfile:
                    diff_fits_file = dcmp_file.replace('cut.dcmp', '').replace('dcmp', '') + 'fits'
                    path_tokens = diff_fits_file.split("/")
                    diff_fits_name = path_tokens.pop()
                    im_path = "/".join(path_tokens)

                    name_tokens = diff_fits_name.split("_")

                    # remove the template portion...
                    name_tokens.pop()
                    name_tokens.pop()
                    name_tokens.pop()

                    im_name = "_".join(name_tokens)
                    im_file = im_path + "/" + im_name + ".sw.fits"

                    csvfile.write('ds9 %s -scale mode 99.5 -region %s %s -scale mode 99.5 -tile -lock frame wcs &\n\n'
                                  % (im_file, unrecovered_reg_file, diff_fits_file))





# MAIN DRIVER SCRIPT
for base_gal_bin_path, gal_bin_subdirs in zip(gal_paths, gal_subdirs):

    print("Processing: %s...\n" % base_gal_bin_path)

    # Spin up objects...
    gal_efficiency_3 = GalEfficiency(base_gal_bin_path, efficiency_magbins, snr=snr, dcmp_type=dcmp_type)
    gal_fakes = [GalFake(base_gal_bin_path, sd) for sd in gal_bin_subdirs]


    # FOR DEBUG PURPOSES
    ds9_commands = "{}/{}".format(gal_efficiency_3.output_dir, "ds9_cmds.txt")
    with open(ds9_commands, 'w') as csvfile:
        csvfile.write("# DS9 Comparison Commands\n")


    # LOAD ALL TEGLON-GENERATED SWOPE TILES
    swope_tile_dict = {}
    swope_tiles = glob.glob("./Fakes/Swope/SwopeTiles/*.txt")
    for st in swope_tiles:
        glade_swope_tile = at.Table.read(st, format='ascii.ecsv')
        model_props = glade_swope_tile.meta['comment']
        dcmp_file = model_props[0].split("=")[1].strip()
        dcmp_base = dcmp_file.replace('.sw.dcmp', '')

        db_id_field_base = st.split("/")[-1].replace('.txt','')
        swope_tile_dict[dcmp_base] = db_id_field_base



    # Iterate over all gal_fakes
    for gf in gal_fakes:
        dcmp_files = os.listdir(gf.dcmp_path)

        # sanity
        for f in dcmp_files:
            if dcmp_type not in f:
                continue

            # region to hold unrecovered stars for this file
            unrecovered_reg_file = f.replace(dcmp_type, "%s.reg" % gf.gal_bin_subdir)
            unrecovered_reg_path = "{}/{}".format(gal_efficiency_3.output_dir, unrecovered_reg_file)

            # # store radial information per fake, per image...
            # radial_reg_file = f.replace(dcmp_type, "%s.radial.reg" % gf.gal_bin_subdir)
            # radial_reg_path = "{}/{}".format(gal_efficiency_3.output_dir, radial_reg_file)


            dcmp_file = "{}/{}".format(gf.dcmp_path, f)
            filename_tokens = f.split(".")
            field_name = filename_tokens[0]
            band = filename_tokens[1]
            utdate = filename_tokens[2].replace("_fake", "")
            img_id = filename_tokens[3][0:4]  # substring the ID out

            img_dcmp_filename = "{}.{}.{}.{}_stch_1".format(field_name, band, utdate, img_id)
            swope_tile_filename_base = swope_tile_dict[img_dcmp_filename]

            # get the corresponding weighted pixels
            pixel_good = at.Table.read("./Fakes/Swope/SwopeDemographics/%s_pixel_good.txt" % swope_tile_filename_base, format='ascii.ecsv')
            weighted_pixel_dict = {}
            for x, y, w in zip(pixel_good['x'], pixel_good['y'], pixel_good['weight']):
                weighted_pixel_dict[(x, y)] = w

            # HACK -- check the band
            if band != current_band:
                continue

            reconstructed_filepath = "{server_img_path}/{field_name}.{band}.{utdate}.{img_id}_stch_1.sw.dcmp".format(
                server_img_path=server_img_path,
                field_name=field_name,
                band=band,
                utdate=utdate,
                img_id=img_id)

            dcmp = txtobj(dcmp_file, cmpheader=True)
            zpt = fits.getval(dcmp_file, 'IZPTMAG')

            fakes = txtobj(gf.fakes_path)
            gf.measure_fakemags(fakes, dcmp, reconstructed_filepath, zpt, gal_efficiency_3.measured_fakes_output,
                                weighted_pixel_dict, unrecovered_reg_file=unrecovered_reg_path, dcmp_file=dcmp_file,
                                unrecovered_range=unrecovered_range, ds9_cmd_file=ds9_commands)










    # FITS
    x0 = 20.0
    a = 0.8
    n = 0.99
    initial_guess = np.asarray([x0, a, n])
    gal_efficiency_3.compute_efficiency(initial_guess)


    # PLOTS
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    ax.plot(model_mags, gal_efficiency_3.computed_lim_mag_model(model_mags),
            'g', label="Galaxy Mag `%s` SNR 3" % base_gal_bin_path)
    ax.errorbar(gal_efficiency_3.central_mag_bins, gal_efficiency_3.weighted_efficiency,
                yerr=[gal_efficiency_3.error_high, gal_efficiency_3.error_low], fmt='kx')

    ax.set_xlabel("mag")
    ax.set_ylabel("Efficiency")
    ax.grid()
    ax.legend()

    output_plot = "{}/{}_{}.png".format(gal_efficiency_3.output_dir, base_gal_bin_path, "SNR_%s" % snr)
    fig.savefig(output_plot, bbox_inches='tight', dpi=300)
    plt.close('all')

    print("Done.")