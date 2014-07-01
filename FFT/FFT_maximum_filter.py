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
sys.path.insert(0, '/home/staff/ggwillc/Desktop/Python_scripts/functions')
import FFT_functions as FFT_functions

startTime = time.time()

# Register driver
#gdal.AllRegister() #<-- useful ont *ly if reading in 
driver = gdal.GetDriverByName('ENVI') ## http://www.gdal.org/formats_list.html
driver.Register()

for file_name in glob("helheim_222a_sample_MAX_ROI_2_butterworth_1_low_pass_fft_output_butterworth_1_filter_50_percent.bin"):
#for file_name in glob("*ELEV*output*butter*1*50_percent.bin"):
#for file_name in glob("*MAX*output*butter*1*50_percent.bin"):

	snip_file_name = file_name.split('.')[0]

	full_snip_1 = file_name.split('_sample')[0]
	print full_snip_1
		
	temp_snip_2 = file_name.split('_pass')[0]
	
	if "MAX" in temp_snip_2:
		full_snip_2 = temp_snip_2.split('MAX_')[1]
		print full_snip_2
	elif "ELEVATION" in temp_snip_2:
		full_snip_2 = temp_snip_2.split('ELEVATION_')[1]
		print full_snip_2
	else:
		"Name snipping error!!"
		os._exit(1)
		
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

	print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	print 'CHECK OUTPUT DIRECTORY EXISTS'
	print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

	#opath = r'/home/staff/ggwillc/Desktop/FFT_image_tests/'
	opath = r'/geog/data/sirius/epsilon/ggwillc/FFT_2D/Helheim/AL_FFT_outputs/small_roi/max_filter/'

	if os.path.isdir(opath):
		print "output_path exists"	
	else:
		print "output_path DOESN'T exist...\n"
		os.makedirs(opath)
		print "...but it does now"		

	print '~~~~~~~~~~~~~~~~~~~~~~~~~~'
	print '~~~~~~~~~~~~~~~~~~~~~~~~~~'
	print 'Processing starts here....'
	print '~~~~~~~~~~~~~~~~~~~~~~~~~~'
	print '~~~~~~~~~~~~~~~~~~~~~~~~~~'

	FFT_surface = image_array

	## Get rolled FFT image object
	plot_title = "Rolled FFT"
	freq = 300
	post = 0.5
	input_x, input_y, magnitude = FFT_functions.magnitude_2D_RETURN(FFT_surface, freq)
	
	def FFT_max_filter_values(freq, post, input_x, input_y, magnitude, kernel = 50):
		print "magnitude.shape"
		print magnitude.shape
		print "magnitude.dtype"
		print magnitude.dtype
		print "magnitude.max"
		print magnitude.max()
		print "magnitude.min"
		print magnitude.min()
		
		print '~~~~~~~~~~~~~~'
		print 'Maximum filter'
		print '~~~~~~~~~~~~~~'

		print "Calculating...."
		method = "_maximum"
		dem_maximum_filter = ndimage.filters.maximum_filter(magnitude,size=(kernel,kernel),mode='reflect')
		filter_max = dem_maximum_filter.max()
		
		'''
		print "magnitude.shape"	
		print dem_maximum_filter.shape
		print "magnitude.dtype"
		print dem_maximum_filter.dtype
		print "Max filter: max"
		
		print filter_max
		print "Max filter: min"
		print dem_maximum_filter.min()
		'''
	
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		print ' Maximum value coordinates '
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		
		b = dem_maximum_filter == (filter_max)
		coords =  np.column_stack(np.where(b))
		#print coords.shape
		max_coord_x = coords[:,0]
		max_coord_y = coords[:,1] 
		
				
		print "############"
		print "     X      " 
		print "############"
		print max_coord_x
		print max_coord_x.shape
		print max_coord_x.dtype

		print "############"
		print "     Y      " 
		print "############"
		print max_coord_y
		print max_coord_y.shape
		print max_coord_y.dtype

		print "##############"
		print " Co-ordinates " 
		print "##############"
		
		
		############
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		print '  Sort max filter surface  '
		print '    Calc. top 10 values 	  '
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		
		'''
		dem_maximum_filter_int = dem_maximum_filter.astype(int)
		indices = dem_maximum_filter_int.ravel().argsort()
		print dem_maximum_filter_int.ravel()[indices[-10:]]
		'''
		
		def get_order_array(a):
			b = np.empty(a.shape, dtype=int)
			for k, row in enumerate(a):
				b[k] = rankdata(-row, method='dense') - 1
			return b
				
		#ordered_array = get_order_array(dem_maximum_filter)
				
		def get_order_array_2(a):
			a_idx = np.argsort(a, axis=-1)[:, ::-1]
			a_sorted = a[np.arange(a.shape[0])[:, None], a_idx]
			a_diff = np.zeros_like(a_sorted, dtype=np.bool)
			a_diff[:, 1:] = a_sorted[:, :-1] != a_sorted[:, 1:]
			a_sorted_ranks = np.cumsum(a_diff, axis=1)
			a_ranks = a_sorted_ranks[np.arange(a_sorted_ranks.shape[0])[:, None], np.argsort(a_idx, axis=1)]
			return a_ranks
		
		ordered_array = get_order_array_2(dem_maximum_filter)	
		print ordered_array
		
		os._exit(1)
		############

		accum_dist_frq = 0.0
		accum_dist_px = 0.0

		# peak of interest (position in maximum array)
		for POI in range(len(max_coord_x)):
		
			x_co = max_coord_x[POI]
			y_co = max_coord_y[POI]
			#mag = magnitude[x_co][y_co] # WRONG : this would get the unfiltered value from the position of the condition required for the max value coords of the filtered surface
			mag = dem_maximum_filter[x_co][y_co]
	
			#print 'Get x and y positions as distances from central origin of [0,0]'
			x_co_frq = max_coord_x - freq/2
			y_co_frq = max_coord_y - freq/2
	
			'''
			print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			print "Calc. diagnonal distance of x,y at position 1 from origin"
			print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			'''
				
			### In frequencies
			x_co_frq_pos = x_co_frq[POI]
			y_co_frq_pos = y_co_frq[POI]
			xy_dist_frq = math.sqrt(x_co_frq_pos**2 + y_co_frq_pos**2)	# relative to the centre of the FFT which is [0,0]
	
			### In pixels
			x_co_px_pos = input_x/x_co_frq_pos
			y_co_px_pos = input_y/y_co_frq_pos
			xy_dist_px = math.sqrt(x_co_px_pos**2 + y_co_px_pos**2)	# relative to the centre of the FFT which is [0,0]
	
			'''
			print "Position: %i" %POI
			print "Magnitude at position:"
			print mag
	
			print "Distance of position:"
			print "%f (frq.)" %(xy_dist_frq)
			print "%f (px)" %(xy_dist_px)
			print "%f (m)" %(xy_dist_px*post)
			'''
			
			accum_dist_frq += xy_dist_frq
			accum_dist_px += xy_dist_px
		
		## Calc. max magnitude mean point distance
			
		mean_dist_frq = accum_dist_frq/len(max_coord_x)
		mean_dist_px = accum_dist_px/len(max_coord_x)
		
		print "Mean spacing of maximum FFT 'spike':"
		print "%f (frq.)" %(mean_dist_frq) 
		print "%f (px)" %(mean_dist_px)
		print "%f (m)" %(mean_dist_px*post)
		
		f = open('FFT_test.txt', 'w')
		f.write("mag,dist.(frq),dist.(px),dist.(m)\n")
		f.write("%f,%f,%f,%f" %(mag,mean_dist_frq,mean_dist_px,mean_dist_px*post))
		f.write("\n")
		f.close()

	FFT_max_filter_values(freq, post, input_x, input_y, magnitude, 50)
		
	#print 'Clear variables' 
	band = None
	image_array = None
	
