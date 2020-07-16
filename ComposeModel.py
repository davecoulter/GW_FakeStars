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


# First step, get all existing curves...
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

class GalEfficiencyModel:
    def __init__(self, m1, m2, x0, a, n):
        self.m1 = m1
        self.m2 = m2
        self.x0 = x0
        self.a = a
        self.n = n

        self.model_x0 = None
        self.model_a = None
        self.model_n = None

    def limiting_mag_model(self, input_mags):
        y = self.n * (1.0 - (erf(self.a * (input_mags - self.x0)) + 1.0) / 2.0)
        return y

    def limiting_mag_model2(self, input_mags):
        y = self.model_n * (1.0 - (erf(self.model_a * (input_mags - self.model_x0)) + 1.0) / 2.0)
        return y

    def get_key(self):
        return "%s - %s" % (self.m1, self.m2)

gal_models = []

for gm in gal_mag_keys:

    efficiencies_path = "./Fakes/Swope/Galaxy_Fakes/%s_%s.cut.dcmp/efficiency_params.txt" % gm

    x0 = None
    a = None
    n = None

    with open(efficiencies_path, 'r') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=' ')
        next(csv_reader) # skip header

        for row in csv_reader:
            x0 = float(row[0])
            a = float(row[1])
            n = float(row[2])

    if x0 is None or a is None or n is None:
        raise Exception("Could not parse efficiency param file: `%s`. Exiting..." % efficiencies_path)

    gal_models.append(GalEfficiencyModel(gm[0], gm[1], x0, a, n))


model_mags = np.linspace(18.0, 25.0, 500)

# Plot all curves
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)

for g in gal_models:
    ax.plot(model_mags, g.limiting_mag_model(model_mags), label="Galaxy Mag `%s` SNR 3" % g.get_key())

ax.yaxis.set_major_locator(MultipleLocator(0.1))
ax.yaxis.set_minor_locator(MultipleLocator(0.02))
ax.xaxis.set_major_locator(MultipleLocator(0.5))
ax.xaxis.set_minor_locator(MultipleLocator(0.1))

ax.set_xlabel("mag")
ax.set_ylabel("Efficiency")
ax.grid(which='both')
ax.legend()

output_plot = "./Fakes/Swope/Galaxy_Fakes/MeasuredGalaxyCurves.png"
fig.savefig(output_plot, bbox_inches='tight', dpi=300)
plt.close('all')

# Zoom in on transition
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)

for g in gal_models:
    ax.plot(model_mags, g.limiting_mag_model(model_mags), label="Galaxy Mag `%s` SNR 3" % g.get_key())

ax.yaxis.set_major_locator(MultipleLocator(0.1))
ax.yaxis.set_minor_locator(MultipleLocator(0.02))
ax.xaxis.set_major_locator(MultipleLocator(0.5))
ax.xaxis.set_minor_locator(MultipleLocator(0.1))

ax.set_xlabel("mag")
ax.set_ylabel("Efficiency")

ax.set_xlim([19.5, 22])
ax.grid(which='both')
ax.legend()

output_plot = "./Fakes/Swope/Galaxy_Fakes/MeasuredGalaxyCurvesZoom.png"
fig.savefig(output_plot, bbox_inches='tight', dpi=300)
plt.close('all')



fig = plt.figure(figsize=(15, 5))
ax_x0 = fig.add_subplot(131)
ax_a = fig.add_subplot(132)
ax_n = fig.add_subplot(133)


x_vals = []
x0_y = []
a_y = []
n_y = []

for g in gal_models:
    mag = (g.m2 - g.m1)/2.0 + g.m1

    ax_x0.plot(mag, g.x0, 'k.')
    ax_a.plot(mag, g.a, 'b.')
    ax_n.plot(mag, g.n, 'r.')

    # hack to remove outlier
    # print(mag)
    # if mag == 14.25 or mag == 14.75 or mag == 15.75:
    #     continue
    # else:
    x_vals.append(mag)
    x0_y.append(g.x0)
    a_y.append(g.a)
    n_y.append(g.n)
    # else:
    #     print("Skip!")

x_vals = np.asarray(x_vals)
print(np.min(x_vals), np.max(x_vals))

x0_m, x0_int = np.polyfit(x_vals, x0_y, 1)
a_m, a_int = np.polyfit(x_vals, a_y, 1)
n_m, n_int = np.polyfit(x_vals, n_y, 1)

x0_f = lambda input_mag: x0_m*input_mag + x0_int
a_f = lambda input_mag: a_m*input_mag + a_int
n_f = lambda input_mag: n_m*input_mag + n_int


ax_x0.plot(x_vals, x0_f(x_vals), 'k--')
ax_a.plot(x_vals, a_f(x_vals), 'b--')
ax_n.plot(x_vals, n_f(x_vals), 'r--')

ax_x0.set_ylabel("x0")
ax_x0.set_xlabel("mag")

ax_a.set_ylabel("a")
ax_a.set_xlabel("mag")

ax_n.set_ylabel("n")
ax_n.set_xlabel("mag")

output_plot = "./Fakes/Swope/Galaxy_Fakes/GalaxyModelParams.png"
fig.savefig(output_plot, bbox_inches='tight', dpi=300)
plt.close('all')


# Plot model
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)

new_mag_model_input = np.linspace(13.25, 17.25, 9)
new_mag_model_output = lambda x0, a, n, input_mags: n * (1.0 - (erf(a * (input_mags - x0)) + 1.0) / 2.0)

for i, m in enumerate(new_mag_model_input):
    x0 = x0_f(m)
    a = a_f(m)
    n = n_f(m)

    gal_models[i].model_x0 = x0
    gal_models[i].model_a = a
    gal_models[i].model_n = n

    y = new_mag_model_output(x0,a,n,model_mags)

    ax.plot(model_mags, y, label="%s" % m)

ax.set_xlabel("mag")
ax.set_ylabel("Efficiency")
ax.grid(which='both')
ax.legend()

output_plot = "./Fakes/Swope/Galaxy_Fakes/GalaxyModel.png"
fig.savefig(output_plot, bbox_inches='tight', dpi=300)
plt.close('all')



# Zoom in on transition
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)

for i, m in enumerate(new_mag_model_input):

    # DEBUG
    # if m != 13.25:
    #     continue

    x0 = x0_f(m)
    a = a_f(m)
    n = n_f(m)

    gal_models[i].model_x0 = x0
    gal_models[i].model_a = a
    gal_models[i].model_n = n

    y = new_mag_model_output(x0,a,n,model_mags)

    ax.plot(model_mags, y, label="%s" % m)

ax.yaxis.set_major_locator(MultipleLocator(0.1))
ax.yaxis.set_minor_locator(MultipleLocator(0.02))
ax.xaxis.set_major_locator(MultipleLocator(0.5))
ax.xaxis.set_minor_locator(MultipleLocator(0.1))

ax.set_xlim([19.5, 22])
ax.set_xlabel("mag")
ax.set_ylabel("Efficiency")
ax.grid(which='both')
ax.legend()

output_plot = "./Fakes/Swope/Galaxy_Fakes/GalaxyModelZoom.png"
fig.savefig(output_plot, bbox_inches='tight', dpi=300)
plt.close('all')





### Test
intial_gm = 15.5
x0 = x0_f(intial_gm)
a = a_f(intial_gm)
n = n_f(intial_gm)
print(x0, a, n)
print(n_m, n_int)
print(a_m, a_int)
print(x0_m, x0_int)






# Plot model residuals
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)

for g in gal_models:

    y1 = g.limiting_mag_model(model_mags)
    y2 = g.limiting_mag_model2(model_mags)

    residuals = y1-y2
    ax.plot(model_mags, residuals, label='%s' % g.get_key())

ax.set_xlabel("mag")
ax.set_ylabel("Residual")
ax.grid(which='both')
ax.legend()

output_plot = "./Fakes/Swope/Galaxy_Fakes/GalaxyModelResidual.png"
fig.savefig(output_plot, bbox_inches='tight', dpi=300)
plt.close('all')


print("Done.")