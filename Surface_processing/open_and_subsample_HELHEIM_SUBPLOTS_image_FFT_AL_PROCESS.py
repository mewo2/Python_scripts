import os
import time as time # for reading in a timer
import numpy as np # maths functions (arrays etc.)
import math
from matplotlib import pyplot as plt # for ploting
from scipy import signal # for convolution function
from scipy import ndimage # for resampling image
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

# 0 = Linux paths | 1 = Windows paths
# Create image representing freq. (1) or space (0)
operating_system = 1 
fft_on = 0

# start timing
startTime = time.time()

variance = []
lag = []

# start timing
#startTime = time.time()
startTime = time.clock()

# Register driver
#gdal.AllRegister() #<-- useful only if reading in 
driver = gdal.GetDriverByName('ENVI') ## http://www.gdal.org/formats_list.html
driver.Register()

# Clear any previous plots
plt.clf()

## For a range of files
#for roi_number in range(1,6):

## For a single file
roi_number = 4
if(roi_number == 4):
	
	# Set file location
	#file_name = r"/geog/data/sirius/epsilon/ggwillc/Maximum_surface_filtering/Helheim/222/HELHEIM_222a_dem_maximum_filter_kernel_239_20_percent_reduction_crevasse_surface_20_percent_reduction"
	#file_name = r"/geog/data/altair/epsilon/ggwillc/AL_ARSF_GRNLND_2013/LiDAR/201a/post_0.5/bin/dem_median_filter_kernel_121_crevasse_surface"
	#file_name = r"/geog/data/sirius/epsilon/ggwillc/Helheim/helheim_lidar_sorting/222a_lidar/bin/222a.helheim_post_0.5m.bin"
	#file_name = r"/geog/data/sirius/epsilon/ggwillc/Helheim/helheim_lidar_sorting/223-_lidar/bin/223-.helheim_post_0.5m.bin"
	#file_name = r'/geog/data/sirius/epsilon/ggwillc/FFT_2D/Helheim/ROI_small/helheim_222a_sample_ELEVATION_ROI_%i.bin' %(roi_number)
	if(operating_system == 1):
		file_name_path = r"C:\Users\ggwillc\Desktop\FFT\FFT_surfaces_binary"
	elif(operating_system == 0):
		#print "Requires a linux path..."
		file_name_path = r''
	
	if(fft_on == 1):
		file_name = "%s\helheim_222a_sample_ELEVATION_ROI_4_no_filtering_fft_output_no_filtering.bin" %(file_name_path)
	elif(fft_on == 0):
		file_name = "%s\helheim_222a_sample_ELEVATION_ROI_4_no_filtering_fft_main_no_filtering.bin" %(file_name_path)
		
	
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

	'''
	print '~~~~~~~~~~~~~~' 
	print 'Subsample 2D array'
	print '~~~~~~~~~~~~~~'

	#image_array_subsample_DATA = image_array[1505:1568, 7079:7145] ## HELHEIM
	image_array_subsample_DATA = image_array[5117:8665, 1209:6028] ## HELHEIM
	'''

	#fig = plt.figure()
	if(fft_on == 1):
		file_name_suffix = "_freq"
		image_array = np.log(image_array)
		plt.imshow(image_array),plt.colorbar()
		title = "ROI %i | cols: %i | rows: %i\n [FREQ IMAGE NOT ROLLED]" %(roi_number, cols, rows)
		plt.title(title)
	elif(fft_on == 0):
		file_name_suffix = "_space"
		image_array = np.log(image_array)
		plt.imshow(image_array),plt.colorbar()
		title = "ROI %i | cols: %i | rows: %i" %(roi_number, cols, rows)
		plt.title(title)
		
	plt.xlabel("Pixel distance (1 px = 0.5m)")
	plt.ylabel("Pixel distance (1 px = 0.5m)")
	
	if(operating_system == 1):
		odir = r'C:\Users\ggwillc\Desktop\FFT\FFT_OUTPUTS_AL_ROUTINE'
	else:
		odir = r'/geog/data/sirius/epsilon/ggwillc/FFT_2D/Helheim/ROI_small'
	
	oput_file = "ROI_%i%s.png" %(roi_number,file_name_suffix)
	output = "%s/%s" %(odir, oput_file) 
	#plt.show()
	plt.savefig(output)
	
	plt.clf()

# End timer
startTime = time.clock()
endTime = time.clock()

print "Program took %.2gs to run" %(endTime - startTime)