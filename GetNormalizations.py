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


working_base_path = "/data/LCO/Swope/workstch/gw190425/1"

images = []
with open("./Fakes/Swope/all_tiles_ascii.txt", 'r') as input_file:
    images += input_file.read().splitlines()

cols = ['img_file', 'num_masked_pix', 'total_pix', 'normalization']
dtype = ['U64', 'i4', 'i4', 'f8']
normalization_table = Table(dtype=dtype, names=cols)

for i, img in enumerate(images):
    t0 = time.time()
    print("... Processing `%s`; [%s/%s]" % (img, i+1, len(images)))
    mask_file = img.replace('fits', 'mask.fits.gz')

    mask_hdu = fits.open("{}/{}".format(working_base_path, mask_file))
    mask_data = mask_hdu[0].data.astype('float')

    num_masked = len(np.where(mask_data == 144.0))
    total_pix = len(mask_data)
    normalization = (1.0 - float(num_masked))/float(total_pix)

    print("\tnum_masked: %s; total pix: %s; norm: %s" % (num_masked, total_pix, normalization))

    normalization_table.add_row([img, num_masked, total_pix, normalization])
    t1 = time.time()
    print("\t... done. %s sec." % (t1-t0))

output_path = "/data/LCO/Swope/logstch/gw190425/1/tile_normalizations.txt"
normalization_table.write(output_path, overwrite=True, format='ascii.ecsv')