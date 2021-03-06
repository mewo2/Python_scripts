'''
#!/usr/bin/python
'''

####
#### MAX filtering
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

site = 'Helheim'
day = '223'

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
print 'Maximum'
print '~~~~~~~~~~~~~~'

print "Calculating...."
method = "_max"
dem_maximum_filter = ndimage.filters.maximum_filter(image_array,size=(kernel,kernel),mode='reflect')
print dem_maximum_filter.shape

filtered_image = dem_maximum_filter
filtered_image_name = "%s_%s_dem%s_filter_kernel_%i" %(site, day, method, kernel)

print '~~~~~~~~~~~~~~'
print 'Surface differencing 10% less'
print '~~~~~~~~~~~~~~'

max_difference_array = image_array - filtered_image
max_difference_array_10_temp = (max_difference_array/100)*10
max_difference_array = max_difference_array - max_difference_array_10_temp

max_difference_name = filtered_image_name + "_max_diff_10_percent_reduction"
print max_difference_array.shape
print max_difference_array.size

'''
print '~~~~~~~~~~~~~~'
print 'Surface differencing 20% less'
print '~~~~~~~~~~~~~~'

max_difference_array = image_array - filtered_image
max_difference_array_20_temp = (max_difference_array/100)*20
max_difference_array_20 = max_difference_array - max_difference_array_20_temp

max_difference_name = filtered_image_name + "_max_diff_20_percent_reduction"
print max_difference_array.shape
print max_difference_array.size
'''
		
print '~~~~~~~~~~~~~~'
print 'Create output directory'
print '~~~~~~~~~~~~~~'
output_directory = r'/geog/data/sirius/epsilon/ggwillc/Maximum_surface_filtering/%s/%s/difference_surfaces/' %(site, day)
create_output_directory(output_directory)

print '~~~~~~~~~~~~~~'
print 'Export the array subtractions as a new ENVI binary file'
print '~~~~~~~~~~~~~~'

print '~~~~~~~~'
print 'max_difference_array output'
print '~~~~~~~~'
max_diff_output_name = "%s%s.bin" %(output_directory, max_difference_name)
max_diff_array_structure = array_to_raster_format(max_diff_output_name,cols_1,rows_1,bands_1,GDT_Float32)
max_diff_array_structure_proj = project_raster(max_diff_array_structure, geotransform_1, projection_1)
max_diff_output = data_to_raster(max_diff_array_structure_proj, max_difference_array)

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
