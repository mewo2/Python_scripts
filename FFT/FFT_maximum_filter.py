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
import FFT_filter_functions as FFT_filter_functions

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
	
	#print magnitude.shape
	pos_x_coord = 400.
	pos_y_coord = 400.
	frq = 300
	angle = FFT_filter_functions.FFT_max_value_bearings_NORTH_FRQ_LESS_180(pos_x_coord,pos_y_coord,frq)
	print "Angle: %f degrees" %angle
	 
	#FFT_filter_functions.FFT_max_filter_values(freq, post, input_x, input_y, magnitude, 50)
		
	#print 'Clear variables' 
	band = None
	image_array = None
	
