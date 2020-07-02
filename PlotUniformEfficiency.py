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

path_formatter = "{}/{}"

server_img_path = "/data/LCO/Swope/workstch/gw190425/1"
base_path = "./Fakes"

measured_uniform_fakes_file = "measuredfakemags.txt"
measured_uniform_fakes_path = path_formatter.format(base_path, measured_uniform_fakes_file)

uniform_efficiencies_file3 = "efficiencies3.txt"
uniform_efficiencies_file10 = "efficiencies10.txt"
uniform_efficiencies_path3 = path_formatter.format(base_path, uniform_efficiencies_file3)
uniform_efficiencies_path10 = path_formatter.format(base_path, uniform_efficiencies_file10)

uniform_1 = "gw190425_fake_1_gw190425tmpl"
uniform_2 = "gw190425_fake_2_gw190425tmpl"
fakes_file = "fakemags.txt"
uniform_1_path = path_formatter.format(base_path, uniform_1)
uniform_2_path = path_formatter.format(base_path, uniform_2)
uniform_1_fakes = path_formatter.format(uniform_1_path, fakes_file)
uniform_2_fakes = path_formatter.format(uniform_2_path, fakes_file)


def measure_fakemags(fake_txtobj, dcmp_txtobj, reconstructed_filepath, zpt, output_file):
    for x, y, mag, mag_gal, mag_gal_GLADE, GLADE_gal_id \
            in zip(fake_txtobj.x[fake_txtobj.dcmpfile == reconstructed_filepath],
                   fake_txtobj.y[fake_txtobj.dcmpfile == reconstructed_filepath],
                   fake_txtobj.mag[fake_txtobj.dcmpfile == reconstructed_filepath],
                   fake_txtobj.mag_gal[fake_txtobj.dcmpfile == reconstructed_filepath],
                   fake_txtobj.mag_gal_GLADE[fake_txtobj.dcmpfile == reconstructed_filepath],
                   fake_txtobj.GLADE_gal_id[fake_txtobj.dcmpfile == reconstructed_filepath]):

        sep = np.sqrt((x - dcmp_txtobj.Xpos) ** 2. + (y - dcmp_txtobj.Ypos) ** 2.)

        # find where fake mag injection is < 2 pixels from DCMP measured source
        if len(np.where(sep < 2)[0]):

            # source confusion: skip if there are multiple matches within 2 pix
            iMin = np.where(sep == np.min(sep))[0]
            if len(iMin) > 1:
                continue

            # Sanity on the values...
            if (dcmp_txtobj.flux[iMin] < 0.0): # reject if flux is < 0
                with open(output_file, 'a') as fout:
                    print('%s %.2f %.2f %.3f %.3f %.3f %.0f %.0f %.0f %0.f' %
                          (reconstructed_filepath, x, y, mag, mag_gal, mag_gal_GLADE, GLADE_gal_id, -99, -99, -99),
                          file=fout)
            else:
                with open(output_file, 'a') as fout:
                    print('%s %.2f %.2f %.3f %.3f %.3f %.0f %.3f %.3f %.3f' %
                          (reconstructed_filepath, x, y, mag, mag_gal, mag_gal_GLADE, GLADE_gal_id,
                           -2.5 * np.log10(dcmp_txtobj.flux[iMin]) + zpt, # calculate mag
                           1.086 * dcmp_txtobj.dflux[iMin] / dcmp_txtobj.flux[iMin], # calculate mag err
                           dcmp_txtobj.flux[iMin] / dcmp_txtobj.dflux[iMin]), # SNR
                          file=fout)
        else:
            with open(output_file, 'a') as fout:
                print('%s %.2f %.2f %.3f %.3f %.3f %.0f %.0f %.0f %0.f' %
                      (reconstructed_filepath, x, y, mag, mag_gal, mag_gal_GLADE, GLADE_gal_id, -99, -99, -99),
                      file=fout)

for file1 in os.listdir(uniform_1_path):

    print(file1)

    # get matching files for run 1 and run 2
    file1_path = path_formatter.format(uniform_1_path, file1)
    file2_path = path_formatter.format(uniform_2_path, file1)

    # skip for now if there aren't two matching files...
    if os.path.exists(file2_path):

        filename_tokens = file1.split(".")
        field_name = filename_tokens[0]
        band = filename_tokens[1]
        utdate = filename_tokens[2].replace("_fake", "")
        img_id = filename_tokens[3][0:4] # substring the ID out

        reconstructed_filepath = "{server_img_path}/{field_name}.{band}.{utdate}.{img_id}_stch_1.sw.dcmp" \
            .format(server_img_path=server_img_path, field_name=field_name, band=band, utdate=utdate, img_id=img_id)

        dcmp_1 = txtobj(file1_path, cmpheader=True)
        zpt_1 = fits.getval(file1_path, 'IZPTMAG')

        dcmp_2 = txtobj(file2_path, cmpheader=True)
        zpt_2 = fits.getval(file2_path, 'IZPTMAG')

        fakes_1 = txtobj(uniform_1_fakes)
        fakes_2 = txtobj(uniform_2_fakes)

        # Initialize output file
        with open(measured_uniform_fakes_path, 'w') as fout:
            print('# dcmpfile x y sim_mag mag_gal mag_gal_GLADE GLADE_gal_id detmag detmagerr snr', file=fout)

        measure_fakemags(fakes_1, dcmp_1, reconstructed_filepath, zpt_1, measured_uniform_fakes_path)
        measure_fakemags(fakes_2, dcmp_2, reconstructed_filepath, zpt_2, measured_uniform_fakes_path)

    break


bin_size = 0.2 # mag
bright = 18
dim = 25
magbins = np.linspace(18, 25, (dim-bright)/bin_size)
measuredfakemags = txtobj(measured_uniform_fakes_path)

with open(uniform_efficiencies_path3, 'w') as fout:
    print('# mag deteff N', file=fout)

    mag_sims = measuredfakemags.sim_mag
    detmags = measuredfakemags.detmag

    for m1, m2 in zip(magbins[:-1], magbins[1:]):

        iDet = np.where((mag_sims[measuredfakemags.snr >= 3.0] >= m1) &
                        (mag_sims[measuredfakemags.snr >= 3.0] <= m2) &
                        (detmags[measuredfakemags.snr >= 3.0] != -99))[0]
        iAll = np.where((mag_sims >= m1) & (mag_sims <= m2))[0]

        if len(iAll):
            print('%.2f %.3f %s' % ((m1 + m2) / 2., len(iDet) / float(len(iAll)), len(iAll)), file=fout)
        else:
            print('%.2f %.3f %s' % ((m1 + m2) / 2., np.nan, len(iAll)), file=fout)

with open(uniform_efficiencies_path10, 'w') as fout:
    print('# mag deteff N', file=fout)

    mag_sims = measuredfakemags.sim_mag
    detmags = measuredfakemags.detmag

    for m1, m2 in zip(magbins[:-1], magbins[1:]):

        iDet = np.where((mag_sims[measuredfakemags.snr >= 10.0] >= m1) &
                        (mag_sims[measuredfakemags.snr >= 10.0] <= m2) &
                        (detmags[measuredfakemags.snr >= 10.0] != -99))[0]
        iAll = np.where((mag_sims >= m1) & (mag_sims <= m2))[0]

        if len(iAll):
            print('%.2f %.3f %s' % ((m1 + m2) / 2., len(iDet) / float(len(iAll)), len(iAll)), file=fout)
        else:
            print('%.2f %.3f %s' % ((m1 + m2) / 2., np.nan, len(iAll)), file=fout)



efficiency3 = []
number3 = []
mag_bins3 = []
error_high3 = []
error_low3 = []

efficiency10 = []
number10 = []
mag_bins10 = []
error_high10 = []
error_low10 = []

with open(uniform_efficiencies_path3, 'r') as fin:
    csv_reader = csv.reader(fin, delimiter=' ')
    next(csv_reader)

    for row in csv_reader:
        mag_bins3.append(float(row[0]))
        efficiency3.append(float(row[1]))
        number3.append(float(row[2]))

for index, e in enumerate(efficiency3):
    yes = e * number3[index]
    no = (1.0 - e) * number3[index]
    ci = binomial(yes, no)

    error_high3.append(ci[0])
    error_low3.append(ci[1])


with open(uniform_efficiencies_path10, 'r') as fin:
    csv_reader = csv.reader(fin, delimiter=' ')
    next(csv_reader)

    for row in csv_reader:
        mag_bins10.append(float(row[0]))
        efficiency10.append(float(row[1]))
        number10.append(float(row[2]))

for index, e in enumerate(efficiency10):
    yes = e * number10[index]
    no = (1.0 - e) * number10[index]
    ci = binomial(yes, no)

    error_high10.append(ci[0])
    error_low10.append(ci[1])

# Limiting mag model uses Eq 12 from:
#   http://www.cfht.hawaii.edu/Science/CFHLS/T0007/T0007-docsu12.html
#   I have removed the 100x prefactor to keep in decimal percent, and introduced a normalization `n` that allows
#   for the maximum efficiency to be different than 1.0
#   Parameters:
#       x0 = turnover in function
#       a = slope at turnover
#       n = normalization

### Version that allows a shift in input mags...
# def limiting_mag_model(input_mags, x0, a, n, d_mag=0.0):
#     delta_mags = np.zeros(len(input_mags)) + d_mag
#     y = n*(1.0 - (erf(a*((input_mags + delta_mags) - x0)) + 1.0)/2.0) # Removed the 100x prefactor
#     return y

def limiting_mag_model(input_mags, x0, a, n):
    y = n*(1.0 - (erf(a*(input_mags - x0)) + 1.0)/2.0)
    return y

def chi_square_s_curve(params, mag_bins=None, obs_efficiency=None, errors=None):
    x0, a, n = params
    model_efficiency = limiting_mag_model(mag_bins, x0, a, n)
    chi_square = np.sum(((model_efficiency - obs_efficiency) / errors) ** 2.0)
    return chi_square

model_mags = np.linspace(18.0, 25.0, 500)

x0=20.0
a=0.8
n=0.99
initial_guess=np.asarray([x0, a, n])

# Estimate errors as average for now...
errors3 = np.mean([error_high3, error_low3], axis=0)
result3 = minimize(chi_square_s_curve, initial_guess, args=(mag_bins3, efficiency3, errors3))
fitted_params3 = result3.x
x0_3_fit = fitted_params3[0]
a_3_fit = fitted_params3[1]
n_3_fit = fitted_params3[2]
fitted_efficiency3 = limiting_mag_model(model_mags, x0_3_fit, a_3_fit, n_3_fit)

errors10 = np.mean([error_high10, error_low10], axis=0)
result10 = minimize(chi_square_s_curve, initial_guess, args=(mag_bins10, efficiency10, errors10))
fitted_params10 = result10.x
x0_10_fit = fitted_params10[0]
a_10_fit = fitted_params10[1]
n_10_fit = fitted_params10[2]
fitted_efficiency10 = limiting_mag_model(model_mags, x0_10_fit, a_10_fit, n_10_fit)


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)

ax.plot(model_mags, fitted_efficiency10,'g', label="SNR 10")
ax.errorbar(mag_bins10, efficiency10, yerr=[error_high10, error_low10], fmt='b+')

ax.plot(model_mags, fitted_efficiency3,'r', label="SNR 3")
ax.errorbar(mag_bins3, efficiency3, yerr=[error_high3, error_low3], fmt='kx')

ax.set_xlabel("mag")
ax.set_ylabel("Efficiency")
ax.legend()

fig.savefig(path_formatter.format(base_path, "fake_mag_test.png"), bbox_inches='tight', dpi=300)
plt.close('all')