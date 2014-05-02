####
#### Surface trend removal
#### @author Chris 29/11/13
####

#
#NB/
#Rows and columns and x origin are a digit out once read 
#in with gdal compared to the dem_par file.....
#

import os
import time as time # for reading in a timer
import numpy as np # maths functions (arrays etc.)
from osgeo import gdal, gdalconst # for reading in raster
from osgeo.gdalconst import * # for reading in raster
from matplotlib import pyplot as plt # for ploting
from scipy import signal # for convolution function
from scipy import ndimage # for resampling image
from matplotlib import cm # colour mapping

# start timing
startTime = time.time()

# Register driver
#gdal.AllRegister() #<-- useful only if reading in 
driver = gdal.GetDriverByName('ENVI') ## http://www.gdal.org/formats_list.html
driver.Register()

## Set file location
#file_name = r"/geog/data/altair/epsilon/ggwillc/AL_ARSF_GRNLND_2013/LiDAR/201a/post_0.5/bin/201a.wholeGlacier.0.5.bin" # the r escapes numbers (and special characters like capital letters) in the pathname
#file_name = r"/geog/data/altair/epsilon/ggwillc/AL_ARSF_GRNLND_2013/LiDAR/203a/post_0.5/bin/203a.wholeGlacier.0.5.bin" # the r escapes numbers (and special characters like capital letters) in the pathname
file_name = r"/geog/data/sirius/epsilon/ggwillc/Helheim/helheim_lidar_sorting/222a_lidar/bin/222a.helheim_post_0.5m.bin"

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
print image_array.shape

print '~~~~~~~~~~~~~~'
print 'Create filter'
print '~~~~~~~~~~~~~~'
kernel = 29# 239 #121
sizex = kernel 
sizey = kernel 

filter1 = np.ones([kernel,kernel])
filter1 = filter1/(kernel*kernel)
print type(filter1)
print filter1.shape

print '~~~~~~~~~~~~~~'
print 'Median'
print '~~~~~~~~~~~~~~'

method = "_median"
dem_median_filter = signal.medfilt2d(image_array,kernel_size=kernel)
print dem_median_filter.shape

filtered_image = dem_median_filter
filtered_image_name = "HELHEIM_222a_dem_median_filter_kernel_%i" % kernel

print '~~~~~~~~~~~~~~'
print 'Output array creation'
print '~~~~~~~~~~~~~~'

median_difference_array = image_array - filtered_image
median_difference_array_name = filtered_image_name + "_median_difference_array"
print median_difference_array.shape
print median_difference_array.size

def negative_anomaly_array(array):
	for i in xrange(len(array)):
		for j in xrange(len(array)):
			if array[i,j] > 0:
				negative_anomaly_array[i,j] = 0
			else:
				negative_anomaly_array[i,j] = negative_anomaly_array[i,j]
	return negative_anomaly_array

def positive_anomaly_array(array):
	for i in xrange(len(array)):
		for j in xrange(len(array)):
			if array[i,j] < 0:
				positive_anomaly_array[i,j] = 0
			else:
				positive_anomaly_array[i,j] = postive_anomaly_array[i,j]
	return negative_anomaly_array
	
negative_anomaly_array  = negative_anomaly_array(median_difference_array)
negative_anomaly_array_name = filtered_image_name + "_negative_anomaly_array"

positive_anomaly_array  = positive_anomaly_array(median_difference_array)
positive_anomaly_array_name = filtered_image_name + "_positive_anomaly_array"
	
print '~~~~~~~~~~~~~~'
print 'Create output directory'
print '~~~~~~~~~~~~~~'
opath = r'/geog/data/sirius/epsilon/ggwillc/Median_surface_filtering/Helheim/'
if os.path.isdir(opath):
	print "output_path exists"	
else:
	print "output_path DOESN'T exist...\n"
	os.makedirs(opath)
	print "...but it does now"

print '~~~~~~~~~~~~~~'
print 'Export the array subtractions as a new ENVI binary file'
print '~~~~~~~~~~~~~~'

#~~~~~~~~
#Median_difference_array
#~~~~~~~~
output_path_median_difference = opath + median_difference_array_name
# Creates a new raster data source
outDs = driver.Create(output_path_median_difference, cols, rows, bands, gdal.GDT_Float32)
# Write metadata
outDs.SetGeoTransform(inds.GetGeoTransform())
outDs.SetProjection(inds.GetProjection())
#Write raster datasets
#for i in range(1):
outBand = outDs.GetRasterBand(1)
outBand.WriteArray(median_difference_array)

# Close raster file
outDs = None
outBand = None

#~~~~~~~~
#negative_anomaly_array
#~~~~~~~~
output_path_negative_anomaly_array = opath + negative_anomaly_array_name
# Creates a new raster data source
outDs = driver.Create(output_path_negative_anomaly_array, cols, rows, bands, gdal.GDT_Float32)
# Write metadata
outDs.SetGeoTransform(inds.GetGeoTransform())
outDs.SetProjection(inds.GetProjection())
#Write raster datasets
#for i in range(1):
outBand = outDs.GetRasterBand(1)
outBand.WriteArray(negative_anomaly_array)

# Close raster file
outDs = None
outBand = None

#~~~~~~~~
#positive_anomaly_array
#~~~~~~~~
output_path_positive_anomaly_array = opath + positive_anomaly_array_name
# Creates a new raster data source
outDs = driver.Create(output_path_positive_anomaly_array, cols, rows, bands, gdal.GDT_Float32)
# Write metadata
outDs.SetGeoTransform(inds.GetGeoTransform())
outDs.SetProjection(inds.GetProjection())
#Write raster datasets
#for i in range(1):
outBand = outDs.GetRasterBand(1)
outBand.WriteArray(output_path_positive_anomaly_array)

# Close raster file
outDs = None
outBand = None

print '~~~~~~~~~~~~~~'
print 'Clear variables'
print '~~~~~~~~~~~~~~'
 
band = None
image_array = None

print '~~~~~~~~~~~~~~'
print '~~~~~~~~~~~~~~'
# end timing
endTime = time.time()
print 'The script took ' + str(endTime - startTime) + ' seconds'

plt.clf() # clears figures
