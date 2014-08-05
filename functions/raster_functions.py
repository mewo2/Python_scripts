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

def ENVI_raster_binary_to_2d_array(file_name):

	# Register driver
	#gdal.AllRegister() #<-- useful only if reading in 

	driver = gdal.GetDriverByName('ENVI') ## http://www.gdal.org/formats_list.html
	driver.Register()

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
	
		# Set pixel offset.....
		print '~~~~~~~~~~~~~~' 
		print 'Convert image to 2D array'
		print '~~~~~~~~~~~~~~'
		band = inds.GetRasterBand(1)
		image_array = band.ReadAsArray(0, 0, cols, rows)
		image_array_name = file_name
		print type(image_array)
		print shape(image_array)
		
		return inds, cols, rows, bands, originX, originY, pixelWidth, pixelHeight, image_array, image_array_name
	