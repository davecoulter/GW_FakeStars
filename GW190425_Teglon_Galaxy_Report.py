# region imports
import os
import sys
import optparse
import csv
import time
import pickle
from collections import OrderedDict
import itertools
import psutil
import shutil
import urllib.request
import requests
from bs4 import BeautifulSoup
from dateutil.parser import parse
import glob
import gc
import json
import pdb
import re
import pytz
from functools import reduce

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.pyplot import cm
from matplotlib.patches import CirclePolygon
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
from scipy.special import erf
from scipy.optimize import minimize, minimize_scalar
from scipy import stats
from scipy.stats import norm
from scipy.integrate import simps, quad, trapz
from scipy.interpolate import interp1d, interp2d

import astropy as aa
from astropy import cosmology
from astropy.coordinates import Distance
from astropy.coordinates.angles import Angle
from astropy.cosmology import z_at_value
from astropy import units as u
import astropy.coordinates as coord
from astropy.table import Table
from astropy.time import Time, TimeMJD

import shapely as ss
from shapely.ops import transform as shapely_transform
from shapely.geometry import Point
from shapely.ops import linemerge, unary_union, polygonize, split
from shapely import geometry

import healpy as hp
from ligo.skymap import distance
# endregion


# Load all tiles, and report:
    # number of unique Swope fields (check that there are 51)
    # number of unique Teglon galaxies within those fields
    # distribution of galaxy mags


dir_formatter = "{}/{}"
base_dir = "SwopeTiles"


swope_tile_files = []
all_B = []
for file_index, file in enumerate(os.listdir(base_dir)):
    file_path = dir_formatter.format(base_dir, file)
    swope_tile_files.append(file_path)

    results_table = Table.read(file_path, format='ascii.ecsv')

    gal_ID = list(results_table['Galaxy_ID'])
    gal_RA = list(results_table['Galaxy_RA'])
    gal_DEC = list(results_table['Galaxy_Dec'])
    PGC = list(results_table['PGC'])
    name_GWGC = list(results_table['Name_GWGC'])
    name_HyperLEDA = list(results_table['Name_HyperLEDA'])
    name_2MASS = list(results_table['Name_2MASS'])
    name_SDSS_DR12 = list(results_table['Name_SDSS_DR12'])
    prob_4D = list(results_table['Galaxy_4D_Prob'])
    z = list(results_table['z'])
    z_dist = list(results_table['z_Dist'])
    z_dist_err = list(results_table['z_Dist_Err'])
    B = list(results_table['B'])
    K = list(results_table['K'])

    all_B += B

fig = plt.figure(figsize=(10, 10), dpi=1000)
ax = fig.add_subplot(111)

ax.hist(all_B)

# ax.set_xlim([])
# ax.set_ylim([])

ax.set_ylabel('Count', fontsize=32, labelpad=9.0)
ax.set_xlabel('mag', fontsize=32)

fig.savefig('galaxy_B_dist.png', bbox_inches='tight')
plt.close('all')
print("... Done.")


