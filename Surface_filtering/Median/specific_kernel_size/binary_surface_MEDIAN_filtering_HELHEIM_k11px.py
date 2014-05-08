'''
#!/usr/bin/python
'''

####
#### MEDIAN filtering
#### @author Chris 07/05/14
####

#
#NB/
#Rows and columns and x origin are a digit out once read 
#in with gdal compared to the dem_par file.....
#

import sys
import os
import time as time # for reading in a timer
import numpy as np # maths functions (arrays etc.)
from osgeo import gdal, gdalconst # for reading in raster
from osgeo.gdalconst import * # for reading in raster
from matplotlib import pyplot as plt # for ploting
from scipy import signal # for convolution function
from scipy import ndimage # for resampling image
from matplotlib import cm # colour mapping

#import home.staff.ggwillc.Desktop.Python_scripts.python_functions.raster_functions

# start timing
startTime = time.time()

# Register driver
#gdal.AllRegister() #<-- useful only if reading in 
driver = gdal.GetDriverByName('ENVI') ## http://www.gdal.org/formats_list.html
driver.Register()

## Set file location
#file_name = r"/geog/data/altair/epsilon/ggwillc/AL_ARSF_GRNLND_2013/LiDAR/201a/post_0.5/bin/201a.wholeGlacier.0.5.bin" # the r escapes numbers (and special characters like capital letters) in the pathname
#file_name = r"/geog/data/altair/epsilon/ggwillc/AL_ARSF_GRNLND_2013/LiDAR/203a/post_0.5/bin/203a.wholeGlacier.0.5.bin" # the r escapes numbers (and special characters like capital letters) in the pathname

#file_name = r"/geog/data/sirius/epsilon/ggwillc/Helheim/helheim_lidar_sorting/222a_lidar/bin/222a.helheim_post_0.5m.bin"
file_name = r"/geog/data/sirius/epsilon/ggwillc/Helheim/helheim_lidar_sorting/223-_lidar/bin/223-.helheim_post_0.5m.bin"

'''
print "Data collection site (e.g. KNS or Helheim etc.):"
sys.stdout.flush()
site = raw_input()

print "Data collection day (julian e.g. 222):"
sys.stdout.flush()
day = raw_input()
'''

site = 'KNS'
day = '201'

#~~~~~~~~~~~~~~~~
#Raster Functions
#~~~~~~~~~~~~~~~~


def open_file_gdal(file_name):
	inds = gdal.Open(file_name, GA_ReadOnly)

	if inds is None:
		print "Really sorry Sir but I couldn't open this blasted file: " + file_name
		print '\nPerhaps you need an ENVI .hdr file? If so, just open the binary up in ENVI and one will be created for you!'
		os._exit(1)
	else:
		print "%s opened successfully" %file_name
	return inds
	

def image_size(inds):
	cols = inds.RasterXSize
	rows = inds.RasterYSize
	bands = inds.RasterCount
	print "columns: %i" %cols
	print "rows: %i" %rows
	print "bands: %i" %bands
	return cols, rows, bands


def georef_info(inds):
	geotransform = inds.GetGeoTransform()
	originX = geotransform[0]
	originY = geotransform[3]
	pixelWidth = geotransform[1]
	pixelHeight = geotransform[5]

	print "origin x: %i" %originX
	print "origin y: %i" %originY
	print "width: %2.2f" %pixelWidth
	print "height: %2.2f" %pixelHeight
	
	projection = inds.GetProjection()
	return geotransform, projection

def convert_to_2darray(inds, cols, rows, file_name):
	band = inds.GetRasterBand(1)
	datatype = band.DataType
	image_array = band.ReadAsArray(0, 0, cols, rows)
	image_array_name = file_name
	print type(image_array)
	print image_array.shape
	return image_array, image_array_name, band, datatype
	
################

def create_output_directory(opath):
	if os.path.isdir(opath):
		print "output_path exists"	
	else:
		print "output_path DOESN'T exist...\n"
		os.makedirs(opath)
		print "...but it does now"
	
################

def array_to_raster_format(array, cols, rows, bands, GDALdatatype):
	raster = driver.Create(array, cols, rows, bands, GDALdatatype)
	return raster


def project_raster(raster_projected, geotransform, projection):
	raster_projected.SetGeoTransform(geotransform)
	raster_projected.SetProjection(projection)		
	geotransform_raster = raster_projected.GetGeoTransform()
	
	originX = geotransform_raster[0]
	originY = geotransform_raster[3]
	pixelWidth = geotransform_raster[1]
	pixelHeight = geotransform_raster[5]
	
	print "origin x: %i" %originX
	print "origin y: %i" %originY
	print "width: %2.2f" %pixelWidth
	print "height: %2.2f" %pixelHeight
	return raster_projected
	
	
def data_to_raster(raster, array):
	raster_data = raster.GetRasterBand(1)
	print "raster_data: %s" %(raster_data)
	raster_data.WriteArray(array, 0, 0)	
	return raster_data
	
	
def negative_anomaly_array(array):
	negative_anomaly_array = array.copy()
	negative_anomaly_array[array > 0] = 0
	return negative_anomaly_array


def positive_anomaly_array(array):
	positive_anomaly_array = array.copy()
	positive_anomaly_array[array < 0] = 0
	return positive_anomaly_array


print '~~~~~~~~~~~~~~'
print 'Open file'
print '~~~~~~~~~~~~~~'

file_1 = open_file_gdal(file_name)

print '~~~~~~~~~~~~~~'
print 'Get image size'
print '~~~~~~~~~~~~~~'

cols_1, rows_1, bands_1 = image_size(file_1)

print '~~~~~~~~~~~~~~'
print 'Get georeference information'
print '~~~~~~~~~~~~~~'

geotransform_1, projection_1 = georef_info(file_1)

print '~~~~~~~~~~~~~~' 
print 'Convert image to 2D array'
print '~~~~~~~~~~~~~~'

image_array, image_array_name, bands, datatype = convert_to_2darray(file_1, cols_1, rows_1, file_name)
print "data_type: %s" %(datatype)

print '~~~~~~~~~~~~~~'
print 'Create filter'
print '~~~~~~~~~~~~~~'
kernel = 11# 239 #121
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
	
negative_anomaly_array  = negative_anomaly_array(median_difference_array)
negative_anomaly_array_name = filtered_image_name + "_negative_anomaly_array"

positive_anomaly_array  = positive_anomaly_array(median_difference_array)
positive_anomaly_array_name = filtered_image_name + "_positive_anomaly_array"
	
print '~~~~~~~~~~~~~~'
print 'Create output directory'
print '~~~~~~~~~~~~~~'
output_directory = r'/geog/data/sirius/epsilon/ggwillc/Median_surface_filtering/%s/%s/difference_surfaces/' %(site, day)
create_output_directory(output_directory)

print '~~~~~~~~~~~~~~'
print 'Export the array subtractions as a new ENVI binary file'
print '~~~~~~~~~~~~~~'


#difference_raster_format = array_to_raster_format('difference_surface.bin', cols_1, rows_1,  bands_1, 'GDT_Float32') # datatypemust be a GDALDataType..

print '~~~~~~~~'
print 'Median_difference_array output'
print '~~~~~~~~'
median_diff_output_name = "%smedian_diff_array_kernel_%i.bin" %(output_directory, kernel)
median_diff_array_structure = array_to_raster_format(median_diff_output_name,cols_1,rows_1,bands_1,GDT_Float32)
median_diff_array_structure_proj = project_raster(median_diff_array_structure, geotransform_1, projection_1)
median_diff_output = data_to_raster(median_diff_array_structure_proj, median_difference_array)

print '~~~~~~~~'
print 'negative_anomaly_array output'
print '~~~~~~~~'
negative_anomaly_output_name = "%snegative_anomaly_array_kernel_%i.bin" %(output_directory, kernel) 
negative_anomaly_array_structure = array_to_raster_format(negative_anomaly_output_name,cols_1,rows_1,bands_1,GDT_Float32)
negative_anomaly_array_structure_proj = project_raster(negative_anomaly_array_structure, geotransform_1, projection_1)
negative_anomaly_output = data_to_raster(negative_anomaly_array_structure_proj, negative_anomaly_array)

print '~~~~~~~~'
print 'positive_anomaly_array output'
print '~~~~~~~~'
positive_anomaly_output_name = "%spositive_anomaly_array_kernel_%i.bin" %(output_directory, kernel) 
positive_anomaly_array_structure = array_to_raster_format(positive_anomaly_output_name,cols_1,rows_1,bands_1,GDT_Float32)
positive_anomaly_array_structure_proj = project_raster(positive_anomaly_array_structure, geotransform_1, projection_1)
positive_anomaly_output = data_to_raster(positive_anomaly_array_structure_proj, positive_anomaly_array)

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
