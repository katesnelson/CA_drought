
from __future__ import division
import numpy as np
import matplotlib as plt
import scipy
from scipy import integrate, interpolate, signal
import scipy.stats as stats
import gdal
from gdalconst import *
gdal.UseExceptions()

#all-time farmland flag construction:
fl = gdal.Open('/data/emily/WF/kate/fl.tif', GA_ReadOnly).ReadAsArray()
fl[fl == 255] = 0
out = np.sum(fl, axis = 0)
out[out > 0] = 1
