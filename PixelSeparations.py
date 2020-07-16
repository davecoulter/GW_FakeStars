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


def dec_2_pix(ra_decimal, dec_decimal, dcmp_textobj):
    alpha = np.radians(ra_decimal - dcmp_textobj.CRVAL1)
    delta = np.radians(dec_decimal)

    eta = (1 + np.cos(alpha) * (np.tan(np.radians(dcmp_textobj.CRVAL2))) ** 2) / (
                np.tan(np.radians(dcmp_textobj.CRVAL2)) + \
                (1 / np.tan(delta))) - np.tan(np.radians(dcmp_textobj.CRVAL2)) * np.cos(alpha)

    top = (1.0 / np.tan(np.radians(dcmp_textobj.CRVAL2)) + np.tan(np.radians(dcmp_textobj.CRVAL2)) * np.cos(alpha))
    bottom = (1.0 + 1.0 / (np.tan(delta) * np.tan(np.radians(dcmp_textobj.CRVAL2))))

    epsilon = (top / bottom) * (np.tan(alpha) / np.tan(delta)) * np.cos(np.radians(dcmp_textobj.CRVAL2))

    x = np.degrees(epsilon)
    y = np.degrees(eta)

    ypix = ((dcmp_textobj.CD2_1 / dcmp_textobj.CD1_1) * (x - dcmp_textobj.CD1_2 * dcmp_textobj.CRPIX2) - \
            (y + dcmp_textobj.CD2_2 * dcmp_textobj.CRPIX2)) / (
                       (dcmp_textobj.CD2_1 / dcmp_textobj.CD1_1) * dcmp_textobj.CD1_2 - dcmp_textobj.CD2_2)
    xpix = (x - dcmp_textobj.CD1_2 * (ypix - dcmp_textobj.CRPIX2)) / dcmp_textobj.CD1_1 + dcmp_textobj.CRPIX1

    return xpix, ypix


def pix_2_dec(x_pixel, y_pixel, dcmp_textobj):
    x = dcmp_textobj.CD1_1 * (x_pixel - dcmp_textobj.CRPIX1) + dcmp_textobj.CD1_2 * (y_pixel - dcmp_textobj.CRPIX2)
    y = dcmp_textobj.CD2_1 * (x_pixel - dcmp_textobj.CRPIX1) + dcmp_textobj.CD2_2 * (y_pixel - dcmp_textobj.CRPIX2)

    epsilon = np.radians(x)
    eta = np.radians(y)

    alpha_tan1 = epsilon / np.cos(np.radians(dcmp_textobj.CRVAL2))
    alpha_tan2 = 1 - eta * np.tan(np.radians(dcmp_textobj.CRVAL2))
    alpha = np.arctan(alpha_tan1 / alpha_tan2)

    ra = np.degrees(alpha) + dcmp_textobj.CRVAL1

    delta_tan1 = eta + np.tan(np.radians(dcmp_textobj.CRVAL2)) * np.cos(alpha)
    delta = np.arctan(delta_tan1 / alpha_tan2)
    dec = np.degrees(delta)

    return ra, dec


fits_name = "s005aae16234.r.ut190428.4151_stch_1.sw.fits"
file_path = "/data/LCO/Swope/workstch/gw190425/1"
test_img_path = '%s/%s' % (file_path, fits_name)

segmap_file_name = fits_name.replace('.fits', '.check2.fits')
sgmap_file_path = '%s/%s' % (file_path, segmap_file_name)

sextable = runsex(test_img_path, segmapname=sgmap_file_path, zpt=fits.getval(test_img_path, 'ZPTMAG'))