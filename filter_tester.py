from __future__ import division

import os
import time as time # for reading in a timer
import numpy as np # maths functions (arrays etc.)
import math
from matplotlib import pyplot as plt # for ploting
from scipy import signal # for convolution function
from scipy import ndimage # for resampling image
from scipy.fftpack import fft2
from matplotlib import cm # colour mapping
import matplotlib as matplotlib
import matplotlib.pyplot as plt
import copy as cp
from numpy import *
import datetime
import random
from scipy import ndimage
from scipy import misc
from scipy.stats import rankdata
from time import gmtime, strftime
import pylab as pl
from glob import glob
from osgeo import gdal, gdalconst # for reading in raster
from osgeo.gdalconst import * # for reading in raster
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import sys
sys.path.insert(0, 'C:/Users/ggwillc/Desktop/Python_scripts/functions')
import Filter_functions
reload(Filter_functions)

a = np.random.rand(5,5)
conservative_filtered_array = Filter_functions.conservative_smooth_efficient(a, 3)