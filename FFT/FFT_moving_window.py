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
from time import gmtime, strftime
import pylab as pl
from osgeo import gdal, gdalconst # for reading in raster
from osgeo.gdalconst import * # for reading in raster

# start timing
startTime = time.time()

variance = []
lag = []

# start timing
startTime = time.time()

# Register driver
#gdal.AllRegister() #<-- useful only if reading in 
driver = gdal.GetDriverByName('ENVI') ## http://www.gdal.org/formats_list.html
driver.Register()

# Clear any previous plots
plt.clf()

# Set file location
#file_name = r"/geog/data/altair/epsilon/ggwillc/AL_ARSF_GRNLND_2013/LiDAR/201a/post_0.5/bin/dem_median_filter_kernel_121_crevasse_surface"
#file_name = r"/geog/data/sirius/epsilon/ggwillc/Maximum_surface_filtering/Helheim/222/HELHEIM_222a_dem_maximum_filter_kernel_239_20_percent_reduction_crevasse_surface_20_percent_reduction" 
file_name = r"/geog/data/arcturus/epsilon/ggwillc/Maximum_surface_filtering/Helheim/222/HELHEIM_222a_dem_maximum_filter_kernel_239_20_percent_reduction_crevasse_surface_20_percent_reduction_SUBSAMPLE"

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

print '~~~~~~~~~~~~~~' 
print 'Subsample 2D array'
print '~~~~~~~~~~~~~~'
#image_array_subsample_DATA = image_array[7119:7219, 7219:7319] ##  [i1:i2, j1:j2]
#image_array_subsample_DATA = image_array[3230:4230, 6117:7117] ## HELHEIM
image_array_subsample_DATA = image_array[5117:8665, 1209:6028] ## HELHEIM
print type(image_array_subsample_DATA)
print shape(image_array_subsample_DATA)

print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
print 'CHECK OUTPUT DIRECTORY EXISTS'
print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

winsize = 512 ## FFT is faster for a window to the power of 2 OR a prime number
stepsize = 1000
#nsample_input = 1000000

#opath = r'/geog/data/sirius/epsilon/ggwillc/vario_outputs/subsample_tests/'
#opath = r'/geog/data/sirius/epsilon/ggwillc/vario_outputs/subsample_tests/winsize_%i_stepsize_%i_random_hits_%i/' %(winsize, stepsize, nsample_input)
opath = r'/geog/data/sirius/epsilon/ggwillc/FFT/tests/'

if os.path.isdir(opath):
	print "output_path exists"	
else:
	print "output_path DOESN'T exist...\n"
	os.makedirs(opath)
	print "...but it does now"

	
def FFT2_processing(image_array):
	print "Calculating FFT..."
	FFT2_output = fft2(image_array)
	return FFT2_output
	
def plot_FFT(FFT2_output, pos_ii, pos_jj, opath):
	plt.clf()
	magnitude = np.absolute(FFT2_output) # gives magnitude component of FFT_output
	magnitude[0,0] = 1
	x, y = magnitude.shape
	magnitude = np.roll(magnitude, x//2, 0) # shifts whole image to middle of axis (x//2)
	magnitude = np.roll(magnitude, y//2, 1)
	maxfreq = 50
	magnitude = magnitude[x//2 - maxfreq: x//2 + maxfreq, y//2 - maxfreq: y//2 + maxfreq]
	fig = plt.figure()
	plt.imshow(np.log(magnitude)),plt.colorbar()
	time_stamp = strftime("%H.%M.%S")
	output_filename = opath + 'window_pos_%i_%i_%s.png' %(pos_ii, pos_jj, time_stamp)
	plt.savefig(output_filename)
	#return FFT_plot

## TEST ARRAY (IGNORE ALL OF THE ABOVE APART FROM THE IMPORTS AND USE THIS IF YOU FANCY)
'''
image_array_subsample_DATA = pl.rand(10,10) #np.ndarray(shape = (30,30)) * random.randrange(10,80,1) # random values between 0 and 80
indNan=np.isnan(image_array_subsample_DATA)
image_array_subsample_DATA[indNan]=10.
'''
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET MOVING WINDOW UP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

if winsize > 1 and winsize%2 != 0:
	print "WINDOW SIZE VALID"
else:
	print "WINDOW SIZE INVALID - MUST BE AN ODD NUMBER AND > 1"

nxwind, nywind = image_array_subsample_DATA.shape

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IMPLEMENT LOOP FOR MOVING WINDOW TO WORK THROUGH
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

imgout = [np.zeros((1 + nxwind // stepsize, 1 + nywind // stepsize))]# for i in xrange(nvars)]
window_iteration=0
for ii in range(0, nxwind, stepsize):
	print "ii: %i" %(ii)
	for jj in range(0, nywind, stepsize):
		
		#print "ii: %i" %(ii)
		#print "jj: %i" %(jj)
		
		# CALCULATE MAX AND MIN RANGES OF ROWS AND COLS THAT CAN BE ACCESSED BY THE WINDOW
		imin=max(0,ii-winsize/2) # gets the maximum of either 0 or ii-winsize/2...
		imax=min(nxwind-1,ii+winsize/2)+1
		jmin=max(0,jj-winsize/2)
		jmax=min(nywind-1,jj+winsize/2)+1
				
		### Cell ii,jj is at centre of array
		datwind=image_array_subsample_DATA[imin:imax,jmin:jmax] # IF CELLS ARE OUTSIDE OF THE ARRAY, THEY TAKE THE VALUES OF THE OTHER END OF THE ARRAY
				
		window_iteration += 1
		
		FFT_ouput = FFT2_processing(datwind)
		plot_FFT(FFT_ouput, ii, jj, opath)

print "Cell hits: %i" %(window_iteration)
