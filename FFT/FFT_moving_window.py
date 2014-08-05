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

# start timing
startTime = time.time()

# Register driver
#gdal.AllRegister() #<-- useful only if reading in 
driver = gdal.GetDriverByName('ENVI') ## http://www.gdal.org/formats_list.html
driver.Register()

# Clear any previous plots
plt.clf()

# Set file location
file_name = r"/geog/data/arcturus/epsilon/ggwillc/Maximum_surface_filtering/Helheim/222/HELHEIM_222a_dem_maximum_filter_kernel_239_20_percent_reduction_crevasse_surface_20_percent_reduction_SUBSAMPLE"

# open file
inds = gdal.Open(file_name, GA_ReadOnly)

if inds is None:
	print "Really sorry Sir but I couldn't open this blasted file: " + file_name
	print '\nPerhaps you need an ENVI .hdr file? If so, just open the binary up in ENVI and one will be created for you!'
	os._exit(1)
else:
	print "%s opened successfully" %file_name
	
inds, cols, rows, bands, originX, originY, pixelWidth, pixelHeight, image_array, image_array_name = raster_functions.ENVI_raster_binary_to_2d_array(file_name)	

print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
print 'CHECK OUTPUT DIRECTORY EXISTS'
print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

kernel_size = 512 ## FFT is faster for a window to the power of 2 OR a prime number
stepsize = 1000

opath = r'/geog/data/sirius/epsilon/ggwillc/FFT2/tests2/'

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
	magnitude[0,0] = 1 # magnitude [0,0] is constant - ignore and make 1
	x, y = magnitude.shape
	magnitude = np.roll(magnitude, x//2, 0) # shifts whole image to middle of axis (x//2)
	magnitude = np.roll(magnitude, y//2, 1)
	maxfreq = 50
	magnitude = magnitude[x//2 - maxfreq: x//2 + maxfreq, y//2 - maxfreq: y//2 + maxfreq]  # cuts down to within 50 frequencies
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

if(kernel_size % 2 != 0 and kernel_size >= 3):
	window = np.ones([kernel_size,kernel_size])
elif(kernel_size % 2 == 0 or kernel_size < 3):
	print "kernel is even - it needs to be odd and at least of a value of 3"
	os._exit(1)

nxwind, nywind = FFT_surface.shape

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IMPLEMENT LOOP FOR MOVING WINDOW TO WORK THROUGH
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#imgout = [np.zeros((1 + nxwind // stepsize, 1 + nywind // stepsize))]# for i in xrange(nvars)] #  filled by each step in the loop (the indiviudal elemets in the loop are affected by window size  - the output image is not)
window_iteration=0
for ii in range(0, nxwind, stepsize):
	print "ii: %i" %(ii)
	for jj in range(0, nywind, stepsize):
		
		#print "ii: %i" %(ii)
		#print "jj: %i" %(jj)
		
		# CALCULATE MAX AND MIN RANGES OF ROWS AND COLS THAT CAN BE ACCESSED BY THE WINDOW
		imin=max(0,i-((kernel_size-1)/2)) # gets the maximum of either 0 or i-kernel_size/2...
		imax=min(nxwind-1,i+((kernel_size-1)/2))+1
		jmin=max(0,j-((kernel_size-1)/2))
		jmax=min(nywind-1,j+((kernel_size-1)/2))+1
			
		calc_moving_window_size_x = imax - imin 
		calc_moving_window_size_y = jmax- jmin
		
		### Cell ii,jj is at centre of array
		datwind=image_array_subsample_DATA[imin:imax,jmin:jmax] # IF CELLS ARE OUTSIDE OF THE ARRAY, THEY TAKE THE VALUES OF THE OTHER END OF THE ARRAY
				
		window_iteration += 1
		
		FFT_ouput = FFT2_processing(datwind)
		plot_FFT(FFT_ouput, ii, jj, opath)

print "Cell hits: %i" %(window_iteration)
