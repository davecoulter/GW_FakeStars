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
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


n_gm = lambda gm, m_n, b_n: m_n*gm + b_n
a_gm = lambda gm, m_a, b_a: m_a*gm + b_a
x0_gm = lambda gm, m_x0, b_x0: m_x0*gm + b_x0
lim_mag = lambda pm, gm, m_n, b_n, m_a, b_a, m_x0, b_x0: \
    n_gm(gm, m_n, b_n)*(1.0 - (erf(a_gm(gm, m_a, b_a)*(pm - x0_gm(gm, m_x0, b_x0))) + 1.0)/2.0)

initial_m_n = 0.00383197475833938
initial_b_n = 0.8911605191259779
initial_m_a = 0.03880483131605417
initial_b_a = 0.8674171983964154
initial_m_x0 = 0.07034637796322604
initial_b_x0 = 19.729687822226563
initial_guess = np.asarray([initial_m_n, initial_b_n, initial_m_a, initial_b_a, initial_m_x0, initial_b_x0])

def chi_square_s_curve(params, galaxy_mags=None, psf_mags=None, obs_efficiency=None, errors=None):
    m_n, b_n, m_a, b_a, m_x0, b_x0 = params
    model_efficiency = lim_mag(psf_mags, galaxy_mags, m_n=m_n, b_n=b_n, m_a=m_a, b_a=b_a, m_x0=m_x0, b_x0=b_x0)
    chi_square = np.sum(((model_efficiency - obs_efficiency) / errors) ** 2.0)
    return chi_square

gal_mag_keys = [
    (13.0, 13.5),
    (13.5, 14.0),
    (14.0, 14.5),
    (14.5, 15.0),
    (15.0, 15.5),
    (15.5, 16.0),
    (16.0, 16.5),
    (16.5, 17.0),
    (17.0, 17.5),
]

central_galmags = []
central_psfmag_bins = []
efficiency = []
number = []
error_high = []
error_low = []

for gm in gal_mag_keys:
    efficiencies_path = "./Fakes/Swope/Galaxy_Fakes/%s_%s.cut.dcmp/efficiencies_3.0.txt" % gm
    central_gal_mag = (gm[1] - gm[0]) / 2.0 + gm[0]

    with open(efficiencies_path, 'r') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=' ')
        next(csv_reader) # skip header

        for row in csv_reader:
            central_psfmag_bins.append(float(row[0]))
            efficiency.append(float(row[1]))
            number.append(central_gal_mag)
            central_galmags.append(central_gal_mag)

for eff, num in zip(efficiency, number):
    yes = eff * num
    no = (1.0 - eff) * num
    ci = binomial(yes, no)

    error_high.append(ci[0])
    error_low.append(ci[1])


mean_errors = np.mean([error_high, error_low], axis=0)




minimize_result = minimize(chi_square_s_curve, initial_guess, args=(np.asarray(central_galmags), central_psfmag_bins,
                                                                    efficiency, mean_errors))



fitted_m_n = minimize_result.x[0]
fitted_b_n = minimize_result.x[1]
fitted_m_a = minimize_result.x[2]
fitted_b_a = minimize_result.x[3]
fitted_m_x0 = minimize_result.x[4]
fitted_b_x0 = minimize_result.x[5]
# DEBUG
print(fitted_m_n, fitted_b_n)
print(fitted_m_a, fitted_b_a)
print(fitted_m_x0, fitted_b_x0)



fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)

model_gal_mags = np.linspace(13.25, 17.25, 9)
psf_mag = np.linspace(18.0, 25.0, len(efficiency))
for gm in model_gal_mags:

    y = lim_mag(psf_mag, gm, fitted_m_n, fitted_b_n, fitted_m_a, fitted_b_a, fitted_m_x0, fitted_b_x0)
    ax.plot(psf_mag, y, label="%s" % gm)

ax.yaxis.set_major_locator(MultipleLocator(0.1))
ax.yaxis.set_minor_locator(MultipleLocator(0.02))
ax.xaxis.set_major_locator(MultipleLocator(0.5))
ax.xaxis.set_minor_locator(MultipleLocator(0.1))

ax.set_xlim([19.5, 22])
ax.set_xlabel("mag")
ax.set_ylabel("Efficiency")
ax.grid(which='both')
ax.legend()

output_plot = "./Fakes/Swope/Galaxy_Fakes/SimultaneousZoom.png"
fig.savefig(output_plot, bbox_inches='tight', dpi=300)
plt.close('all')