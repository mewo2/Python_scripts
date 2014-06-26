from __future__ import division

import sys
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
from time import gmtime, strftime
import pylab as pl
from glob import glob

from osgeo import gdal, gdalconst # for reading in raster
from osgeo.gdalconst import * # for reading in raster

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import sys
sys.path.insert(0, '/home/staff/ggwillc/Desktop/Python_scripts/functions')
import FFT_functions as FFT_functions

# start timing
startTime = time.time()

# Register driver
#gdal.AllRegister() #<-- useful only if reading in 
driver = gdal.GetDriverByName('ENVI') ## http://www.gdal.org/formats_list.html
driver.Register()

###################################
###################################
## Start loop
###################################
###################################

# Clear any previous plots
plt.clf()

## Set file location
#file_name = r"/geog/data/sirius/epsilon/ggwillc/FFT_2D/Helheim/AL_FFT_outputs/low_pass/heavy_crevasse_ROI/heavy_crevasse_fft_output_tophat.bin"
#snip_file_name = file_name.split('heavy_crevasse_ROI/')[1]

# IF using this, don't forget to tab everything below!
#for file_name in glob("*.bin"):

def usage():
	print 'Usage: '+sys.argv[0]+' <file.bin> <Frequency_resolution>'
	
args_length = len(sys.argv)
print "args_length: %i" %(args_length)

if (args_length >= 3 ):
	print "args length = %f" %args_length
else:
	usage()
	sys.exit()
	
#file_name = "heavy_crevasse_fft_output_tophat.bin"

#file_name = "lesser_crevasse_2_fft_output_gauss_filter_25_percent.bin"
file_name =  str(sys.argv[1])

#frq = 1000

frq =  int(sys.argv[2])
	
snip_file_name = file_name.split('.')[0]
		
# open file
inds = gdal.Open(file_name, GA_ReadOnly)

if inds is None:
	print "Really sorry Sir but I couldn't open this blasted file: " + file_name
	print '\nPerhaps you need an ENVI .hdr file? If so, just open the binary up in ENVI and one will be created for you!'
	os._exit(1)
else:
	print "%s opened successfully" %file_name
		
print '~~~~~~~~~~~~~~'
print 'Get image size'
print '~~~~~~~~~~~~~~'
cols = inds.RasterXSize
rows = inds.RasterYSize
bands = inds.RasterCount

print "columns: %i" %cols
print "rows: %i" %rows
print "bands: %i" %bands

print '~~~~~~~~~~~~~~'
print 'Get georeference information'
print '~~~~~~~~~~~~~~'
geotransform = inds.GetGeoTransform()
originX = geotransform[0]
originY = geotransform[3]
pixelWidth = geotransform[1]
pixelHeight = geotransform[5]

print "origin x: %i" %originX
print "origin y: %i" %originY
print "width: %2.2f" %pixelWidth
print "height: %2.2f" %pixelHeight

print '~~~~~~~~~~~~~~' 
print 'Convert image to 2D array'
print '~~~~~~~~~~~~~~'

band = inds.GetRasterBand(1)
image_array = band.ReadAsArray(0, 0, cols, rows)
image_array_name = file_name
print type(image_array)
print shape(image_array)

FFT_surface = image_array

print '~~~~~~~~~~~~~~' 
print 'Subsample 2D array'
print '~~~~~~~~~~~~~~'
#image_array_subsample_DATA = image_array[7119:7219, 7219:7319] ##  [i1:i2, j1:j2]
#image_array_subsample_DATA = image_array[3230:4230, 6117:7117] ## HELHEIM
#image_array_subsample_DATA = image_array[5117:8665, 1209:6028] ## HELHEIM
#print type(image_array_subsample_DATA)
#print shape(image_array_subsample_DATA)

print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
print 'CHECK OUTPUT DIRECTORY EXISTS'
print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

opath = r'/home/staff/ggwillc/Desktop/FFT_image_tests/'

if os.path.isdir(opath):
	print "output_path exists"	
else:
	print "output_path DOESN'T exist...\n"
	os.makedirs(opath)
	print "...but it does now"
	
FFT_functions.plot_FFT_2D_interactive_z_AXIS_FREQ(FFT_surface, opath, snip_file_name, frq)
#FFT_functions.plot_FFT_2D_interactive_z_AXIS_WAVELENGTHS(FFT_surface, opath, snip_file_name, frq)
#FFT_functions.plot_FFT_2D_interactive_z(FFT_surface, opath)

#plot_FFT_3D(FFT_surface, opath)

print '~~~~~~~~~~~~~~'
print 'Clear variables'
print '~~~~~~~~~~~~~~'

band = None
image_array = None
